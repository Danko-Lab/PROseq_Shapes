# -*- coding: utf-8 -*-

# To do
# Adjust sample size depending on label (percentages would be good) - Done (05/31/20)
# Increase number of steps to 1000 - Done (05/25/20)
# Use feature midpoint instead of start - Done (05/24/20)
# Add version number to each file (could be version directory instead) - Done (05/24/20)

# Version history:

# Label version 1 (LAB_V1) - original version:
# - Using directory depth[0]-window[100]-step[50]-andor[True]
# - 20Kb regions added to each end of annotation
# - Threshold of 1 TPM for filtering active annotations
# - Single label per window

# Label version 2 (LAB_V2) - revisions after defense:
# - Using directory depth[0]-window[100]-step[50]-andor[True]
# - 20Kb regions added to each end of annotation
# - Multiple labels per window
# - High confidence GENCODE annotations (only those that intersect w/ GRO-cap and poly-A sites)

# CNN version 1 - original version:
# - Window size: 8192 bins, bin size: 16bp, batch size: 128

# CNN version 2 - select center of window at random from within feature (instead of always midpoint)
# - Window size: 8192 bins, bin size: 16bp, batch size: 128

# CNN version 3 - 51K windows
# - Window size: 1024 bins, bin size: 50bp, batch size: 128

# CNN version 4 - 51K windows
# - Window size: 1024 bins, bin size: 50bp, batch size: 128
# - Reduce number of hidden layers to 5 (to compensate for smaller input size)

"""discriminator-multiclass_prm.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/11z7yq_283IZRJ3d3ivZp_ZSPr15s3o7M
"""

"""
Do all the Colab stuff up front
Getting access to the google drive in colab for the dataset. 
"""
#from google.colab import drive
#drive.mount('/gdrive')

# Set up globals
import set_up_globals
labelVersion = set_up_globals.labelVersion
CNNVersion = set_up_globals.CNNVersion
data_folder = set_up_globals.data_folder
model_folder = set_up_globals.model_folder
gpuNumber = set_up_globals.gpuNumber

#!pip install pyBigWig
#!pip install pybedtools

#!ls "{data_folder}"
#!ls

import pyBigWig
import pybedtools
import pandas as pd
import numpy as np
import random
random.seed(2)

PREVIOUS_EPOCH = 0
numberOfStepsPerEpoch = 50
numberOfValidationSteps = 10
pickleFile = 'genebody-dataframe.pkl'

LABELS = ['plus-neg', 'minus-neg', 'plus-genebody', 'minus-genebody',
          'plus-aftergene', 'minus-aftergene', 'plus-genestart',
          'minus-genestart', 'plus-geneend', 'minus-geneend', 'tss',
          'plus-stable', 'minus-stable', 'plus-unstable', 'minus-unstable']

totalSamples = 4000000
plus_neg_sample_percent = 0.2
minus_neg_sample_percent = 0.2
plus_genebody_sample_percent = 0.15
minus_genebody_sample_percent = 0.15
plus_aftergene_sample_percent = 0.025
minus_aftergene_sample_percent = 0.025
plus_genestart_sample_percent = 0.025
minus_genestart_sample_percent = 0.025
plus_geneend_sample_percent = 0.025
minus_geneend_sample_percent = 0.025
tss_sample_percent = 0.05
plus_stable_sample_percent = 0.025
minus_stable_sample_percent = 0.025
plus_unstable_sample_percent = 0.025
minus_unstable_sample_percent = 0.025

LABEL_WEIGHTS = {'plus-neg': int(plus_neg_sample_percent  *  totalSamples), 
          'minus-neg': int(minus_neg_sample_percent  *  totalSamples), 
          'plus-genebody': int(plus_genebody_sample_percent  *  totalSamples), 
          'minus-genebody': int(minus_genebody_sample_percent  *  totalSamples),
          'plus-aftergene': int(plus_aftergene_sample_percent  *  totalSamples), 
          'minus-aftergene': int(minus_aftergene_sample_percent  *  totalSamples), 
          'plus-genestart': int(plus_genestart_sample_percent  *  totalSamples),
          'minus-genestart': int(minus_genestart_sample_percent  *  totalSamples), 
          'plus-geneend': int(plus_geneend_sample_percent  *  totalSamples), 
          'minus-geneend': int(minus_geneend_sample_percent  *  totalSamples), 
          'tss': int(tss_sample_percent  *  totalSamples), 
          'plus-stable': int(plus_stable_sample_percent  *  totalSamples), 
          'minus-stable': int(minus_stable_sample_percent  *  totalSamples), 
          'plus-unstable': int(plus_unstable_sample_percent  *  totalSamples), 
          'minus-unstable': int(minus_unstable_sample_percent  *  totalSamples) }


class ProseqFeatureFactory:
    """Data factory class for a specific train/val/test set
    
    Allows for calling as an iterator to obtain the next training point in the
    dataset. This class does feature extraction on the bigwigs, but for now only
    returns a 100,000 dimensional numpy array as the features. Numpy arrays are padded
    with 0's if it bleeds over the edge of the chromosome.
    """
    def __init__(self, plus_bws, minus_bws, dataset, binsize=None, window=2048):
        self._dataset = dataset
        self._plus_bws = plus_bws
        self._minus_bws = minus_bws
        self._binsize = binsize
        self._window = window
        self._ptr = 0
        
    def __iter__(self):
        return self
    
    def __next__(self):
        # If we finish going through the dataset, shuffle and restart from beginning
        if self._ptr >= len(self._dataset):
            raise StopIteration
        row = self._dataset.iloc[self._ptr, :]
        plus_bw = self._plus_bws[row.dataset]
        minus_bw = self._minus_bws[row.dataset]
        
        # Start pulling features from the bigwigs
        # (center on feature mid-point, not the beginning)
        halfwindow_binned = self._window // 2 * self._binsize
        fullwindow_binned = self._window * self._binsize
        # Using feature mid-point instead of start - CNN version 1
        # Using random point within feature, instead of midpoint - CNN version 2
        #midpoint = row.start + ((row.end - row.start) // 2)
        midpoint = row.start + random.randint(0, (row.end - row.start))
        start = max(0, midpoint - halfwindow_binned)
        total_loci = plus_bw.chroms(row.chrom)
        end = min(total_loci, midpoint + halfwindow_binned)
        try:
            plus_arr = plus_bw.values(row.chrom, start, end, numpy=True)
            minus_arr = minus_bw.values(row.chrom, start, end, numpy=True)
        except:
            self._ptr += 1
            return self.__next__()
        # Pad the features if necessary to have feature vectors of length fullwindow_binned
        if len(plus_arr) != fullwindow_binned and midpoint - halfwindow_binned < 0:
            plus_arr = np.pad(plus_arr, ((halfwindow_binned - midpoint), 0), 'constant')
            minus_arr = np.pad(minus_arr, ((halfwindow_binned - midpoint), 0), 'constant')
        if len(plus_arr) != fullwindow_binned and midpoint + halfwindow_binned > total_loci:
            plus_arr = np.pad(plus_arr, (0, midpoint + halfwindow_binned - total_loci), 'constant')
            minus_arr = np.pad(minus_arr, (0, midpoint + halfwindow_binned - total_loci), 'constant')
        # Stack the features so that the first row is the positive reads and second row is neg reads
        data = np.nan_to_num(np.vstack((plus_arr, np.abs(minus_arr))))
        # Scale the data to range between 0 and 1
        if np.max(data) > 0.:
            data = data / np.max(data)
        # Bin the features if necessary
        if self._binsize:
            data = np.add.reduceat(data, np.arange(0, len(data[0]), self._binsize), axis=1)
        
        # Construct one-hot encoded labels
        lbl = np.array([row[lbl] for lbl in LABELS])
        
        self._ptr += 1
        return data, lbl, (row.dataset, row.chrom, row.start, row.end)
    
class ProseqDataFactory:
    """Data factory class for proseq data. Must be called through a with: statement to instantiate file handlers!
    """
    def __init__(self, dataframes, plus_bws, minus_bws, datasets, window=2048, binsize=None, weight_by_col='coverage100000'):
        """Initializes the data factory
        
            Requires:
            dataframes: {'dataset_name': pd.Dataframe} for each dataset. Dataframes should have the columns
                ['chrom', 'start', 'end', 'coverage100', 'coverage1000', 'coverage10000', 'coverage100000', 'label']
            plus_bws/minus_bws: plus and minus maps of {'dataset_name': 'bw_location'} for the datasets
            truth_set: dataset_name of the set to use as ground truth for the autoencoder
            n_positive: number of positive informative sites to pull for each dataset
            n_negative: number of negative informative sites to pull for each dataset
            binsize: the size of bins to aggregate counts over
            
            Example Usage:
            DATA_LOCATIONS = {
                'G1': ('../../data/seq/G1/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G1/G1_plus.bw', '../../data/seq/G1/G1_minus.bw'),
                'G2': ('../../data/seq/G2/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G2/G2_plus.bw', '../../data/seq/G2/G2_minus.bw'),
                'G3': ('../../data/seq/G3/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G3/G3_plus.bw', '../../data/seq/G3/G3_minus.bw'),
                'G4': ('../../data/seq/G4/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G4/G4_plus.bw', '../../data/seq/G4/G4_minus.bw'),
                'G5': ('../../data/seq/G5/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G5/G5_plus.bw', '../../data/seq/G5/G5_minus.bw'),
                'G6': ('../../data/seq/G6/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G6/G6_plus.bw', '../../data/seq/G6/G6_minus.bw'),
                'G7': ('../../data/seq/G7/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/G7/G7_plus.bw', '../../data/seq/G7/G7_minus.bw'),
                'GM': ('../../data/seq/GM/informative_positions/depth[0]-window[100]-step[50]-andor[True]/genebody-dataframe.pkl', '../../data/seq/GM/GM_plus.bw', '../../data/seq/GM/GM_minus.bw'),
            }

            DATAFRAMES = {dset: pd.read_pickle(DATA_LOCATIONS[dset][0]) for dset in DATA_LOCATIONS}
            PLUS_BWS = {dset: DATA_LOCATIONS[dset][1] for dset in DATA_LOCATIONS}
            MINUS_BWS = {dset: DATA_LOCATIONS[dset][2] for dset in DATA_LOCATIONS}

            with ProseqDataFactory(DATAFRAMES, PLUS_BWS, MINUS_BWS, 'G1', binsize=10) as df:
                for arr, label in df.train_factory():
                    # Do something
                    break
        """
        if set(dataframes.keys()) != set(plus_bws.keys()) or set(plus_bws.keys()) != set(minus_bws.keys()):
            raise ValueError('Incompatible dataset names across dataframes, plus_bws and minus_bws')
        self._window = window
        self._binsize = binsize
        self._dataframes = dataframes
        self._plus_bws = {dset: pyBigWig.open(plus_bws[dset]) for dset in plus_bws}
        self._minus_bws = {dset: pyBigWig.open(minus_bws[dset]) for dset in minus_bws}
        self._datasets = datasets
        self._weight_by_col = weight_by_col
        
        for dset in self._dataframes:
            self._dataframes[dset]['dataset'] = dset
        
        self.factory = None
        self.reset_factory()
        
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            return next(self.factory)
        except StopIteration:
            self.reset_factory()
            return next(self.factory)
        
    def reset_factory(self):
        pd_dset = pd.concat(
            [self._select_from_dset(dset, self._weight_by_col)
             for dset in self._dataframes if dset in self._datasets],
            axis=0,
        ).sample(frac=1)
        self.factory = ProseqFeatureFactory(
            self._plus_bws,
            self._minus_bws,
            pd_dset,
            binsize=self._binsize,
            window=self._window
        )
        
    def __enter__(self):
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        for dset in self._plus_bws:
            self._plus_bws[dset].close()
        for dset in self._minus_bws:
            self._minus_bws[dset].close()
    
    def _select_from_dset(self, dataset_name, weight_by_col):
        """Helper to select appropriate datapoints from the dataset
        """
        df = self._dataframes[dataset_name]
        dfs = [df[df[lbl] == 1].sample(n=LABEL_WEIGHTS[lbl], axis=0, replace=True, weights=weight_by_col) for lbl in LABELS]
        
        return pd.concat(dfs, axis=0).sample(frac=1)

from sklearn.preprocessing import normalize

DATA_LOCATIONS = {
    'G1': (data_folder + 'seq/G1/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G1/G1_plus.bw', data_folder + 'seq/G1/G1_minus.bw'),
    'G2': (data_folder + 'seq/G2/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G2/G2_plus.bw', data_folder + 'seq/G2/G2_minus.bw'),
    'G3': (data_folder + 'seq/G3/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G3/G3_plus.bw', data_folder + 'seq/G3/G3_minus.bw'),
    #'G4': (data_folder + 'seq/G4/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G4/G4_plus.bw', data_folder + 'seq/G4/G4_minus.bw'),
    'G5': (data_folder + 'seq/G5/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G5/G5_plus.bw', data_folder + 'seq/G5/G5_minus.bw'),
    'G6': (data_folder + 'seq/G6/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G6/G6_plus.bw', data_folder + 'seq/G6/G6_minus.bw'),
    'G7': (data_folder + 'seq/G7/' + labelVersion + '/' + pickleFile, data_folder + 'seq/G7/G7_plus.bw', data_folder + 'seq/G7/G7_minus.bw'),
    #'GM': (data_folder + 'seq/GM/' + labelVersion + '/' + pickleFile, data_folder + 'seq/GM/GM_plus.bw', data_folder + 'seq/GM/GM_minus.bw'),
}

DATAFRAMES = {dset: pd.read_pickle(DATA_LOCATIONS[dset][0]) for dset in DATA_LOCATIONS}
PLUS_BWS = {dset: DATA_LOCATIONS[dset][1] for dset in DATA_LOCATIONS}
MINUS_BWS = {dset: DATA_LOCATIONS[dset][2] for dset in DATA_LOCATIONS}

BATCH_SIZE = set_up_globals.BATCH_SIZE
WINDOW = set_up_globals.WINDOW
BINSIZE = set_up_globals.BINSIZE
#BATCH_SIZE = 128
#WINDOW = 8192
#BINSIZE = 16

# Note: removed G4 from training samples - does not correlate well with other datasets
def train_data_generator():
    with ProseqDataFactory(DATAFRAMES, PLUS_BWS, MINUS_BWS, ['G3', 'G5', 'G6', 'G7'], binsize=BINSIZE, window=WINDOW, weight_by_col='coverage100000') as df:
        while True:
            features = []
            labels = []
            for _ in range(BATCH_SIZE):
                try:
                    arr, label, loc = next(df)
                except StopIteration:
                    return
                features.append(arr)
                labels.append(label)
            features = np.expand_dims(np.swapaxes(np.swapaxes(np.dstack(features), 0, 2), 1, 2), axis=3)
            labels = np.vstack(labels)
            yield features, labels
            
def val_data_generator():
    with ProseqDataFactory(DATAFRAMES, PLUS_BWS, MINUS_BWS, ['G2'], binsize=BINSIZE, window=WINDOW, weight_by_col=None) as df:
        while True:
            features = []
            labels = []
            for _ in range(BATCH_SIZE):
                try:
                    arr, label, loc = next(df)
                except StopIteration:
                    return
                features.append(arr)
                labels.append(label)
            features = np.expand_dims(np.swapaxes(np.swapaxes(np.dstack(features), 0, 2), 1, 2), axis=3)
            labels =  np.vstack(labels)
            yield features, labels

from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Activation, Conv2D, ZeroPadding2D, MaxPooling2D, Flatten, Activation, Dropout, Reshape, BatchNormalization
from tensorflow.keras.callbacks import TensorBoard, ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K

# Force keras to only use first GPU
import os
os.environ["CUDA_VISIBLE_DEVICES"]=gpuNumber

OUTPUT_FOLDER = data_folder + 'models/' + labelVersion + '_' + CNNVersion + '/' + model_folder
# Don't forget to change this when setting the directory for tensorboard logs

#PREVIOUS_EPOCH = 0  # Now setting at the top
LOAD_MODEL_FROM = OUTPUT_FOLDER + '/weights-{}.hdf5'.format(str(PREVIOUS_EPOCH).zfill(4)) if PREVIOUS_EPOCH else ''

DROPOUT = 0.2

if LOAD_MODEL_FROM:
    model = load_model(LOAD_MODEL_FROM)
else:
    model = Sequential()
    input_shape = (2, WINDOW, 1)

#     for i in range(5):
#         model.add(Conv2D(32, kernel_size=(2, 4),
#                          strides=(1, 1), padding='same', input_shape=input_shape))
#         model.add(BatchNormalization())
#         model.add(Activation('relu'))
#         if DROPOUT != 0:
#             model.add(Dropout(DROPOUT))

#     model.add(Conv2D(32, kernel_size=(2, 4), strides=(1, 1), padding='same'))
#     model.add(BatchNormalization())
#     model.add(Activation('relu'))
#     model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
#     if DROPOUT != 0:
#         model.add(Dropout(DROPOUT))

    model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same', input_shape=input_shape))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))
    
    model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))
    
    model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))
        
    model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))
        
    #model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    #model.add(BatchNormalization())
    #model.add(Activation('relu'))
    #model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    #if DROPOUT != 0:
    #    model.add(Dropout(DROPOUT))
        
    #model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    #model.add(BatchNormalization())
    #model.add(Activation('relu'))
    #model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    #if DROPOUT != 0:
    #    model.add(Dropout(DROPOUT))
        
    #model.add(Conv2D(32, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    #model.add(BatchNormalization())
    #model.add(Activation('relu'))
    #model.add(MaxPooling2D(pool_size=(1, 2), strides=(1, 2)))
    #if DROPOUT != 0:
    #    model.add(Dropout(DROPOUT))
    
    model.add(Conv2D(16, kernel_size=(3, 8), strides=(1, 1), padding='same'))
    model.add(BatchNormalization())
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2), strides=(1, 2)))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))

    model.add(Flatten())
    model.add(Dense(256, activation='relu'))
    if DROPOUT != 0:
        model.add(Dropout(DROPOUT))
    model.add(Dense(15, activation='sigmoid'))

    opt = Adam(lr=0.0001)

    model.compile(optimizer=opt, loss='binary_crossentropy', metrics=['accuracy'])

#print(model.summary())

# Clear any logs from previous runs
#!rm -rf "{OUTPUT_FOLDER}/logs/"

# Commented out IPython magic to ensure Python compatibility.
#!kill 2202

# Load the TensorBoard notebook extension
# %load_ext tensorboard

from tensorflow.keras.callbacks import TensorBoard, ModelCheckpoint
import datetime

finished_training = False
save_freq = 10

logs_common = '/logs/fit/'
logsfit = '"' + OUTPUT_FOLDER + logs_common + '"'
log_dir = OUTPUT_FOLDER + logs_common + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard = TensorBoard(log_dir=log_dir, histogram_freq=1, write_graph=True, write_images=True)

# Include the epoch in the file name (uses `str.format`)
#checkpoint_path = OUTPUT_FOLDER + "/cp-{epoch:04d}.ckpt"
checkpoint_path = OUTPUT_FOLDER + "/weights-{epoch:04d}.hdf5"
checkpoint_dir = os.path.dirname(checkpoint_path)

# Create a callback that saves the model every epoch
checkpointer = ModelCheckpoint(
    filepath=checkpoint_path, 
    verbose=1, 
    save_weights_only=False,
    save_freq='epoch')

##checkpointer = ModelCheckpoint(filepath=OUTPUT_FOLDER + '/weights.hdf5', verbose=1, save_freq=save_freq, save_best_only=True)
# period argument depreciated
#checkpointer = ModelCheckpoint(filepath=OUTPUT_FOLDER + '/weights.hdf5', verbose=1, period=save_freq, save_best_only=True)

# Load the previously saved weights
#latest = tf.train.latest_checkpoint(checkpoint_dir)
#model.load_weights(latest)

# Save the weights using the `checkpoint_path` format
# model.save_weights(checkpoint_path.format(epoch=0))

# Launching the tensorboard.
# !kill 694
# %tensorboard --logdir $logsfit
#%tensorboard --logdir '/gdrive/My Drive/Colab Notebooks/data/multiclass-large-windows/logs/fit/'

# View open TensorBoard instances.
#from tensorboard import notebook
#notebook.list()

h = model.fit(
    train_data_generator(),
    callbacks=[tensorboard, checkpointer],
    steps_per_epoch=numberOfStepsPerEpoch,
    validation_steps=numberOfValidationSteps,
    epochs=4800,
    validation_data=val_data_generator(),
    initial_epoch=PREVIOUS_EPOCH,
)

'''
h = model.fit_generator(
    train_data_generator(),
    callbacks=[tensorboard, checkpointer],
    steps_per_epoch=50,
    validation_steps=10,
    epochs=4800,
    validation_data=val_data_generator(),
    initial_epoch=PREVIOUS_EPOCH,
)
'''

finished_training = True

#checkpoint_path

# Cell to clear GPU memory after training
if finished_training:
    K.clear_session()
    print('Done')
    exit(0)


