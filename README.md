# PROseq_Shapes
Comparison of two neural network methods for transcription unit annotation using precision nuclear run-on (PRO-seq)

Understanding the functional elements of the genome that regulate gene expression provides critical insight into the processes governing cell development and differentiation, and related disease mechanisms. An effective strategy for annotating functional elements is to measure de novo transcription along the genome. However, the relationship between transcription and functional element annotations is far from straightforward, and features within the transcriptional signal are often subtle and challenging to model. Assays such as precision nuclear run-on and sequencing (PRO-seq) are a rich source of transcriptional information and numerous annotation methods have been developed to mine these data. We hypothesize that there is additional information in PRO-seq data that has not been leveraged by previous methods. Here we investigate two computational approaches that attempt to identify and model this information and use it to improve genome annotations across a variety of cell types. Each method employs a neural network since these are well suited to the non-linear nature of the task. The first is a convolutional neural network (CNN) combined with a hidden Markov model (HMM), and the second is a conditional generative adversarial network (cGAN). Comparisons with existing methods that used the same input data showed a performance improvement for the first of our methods (CNN-HMM). The CNN used by this method was able to reliably pick out multiple features from the PRO-seq data, thus providing evidence for our original hypothesis. The cGAN, in contrast, shows great potential for producing fully labeled, simulated PRO-seq data.

## Citation

If you use this code or the resulting assemblies, please cite the following paper:

*Comparison of two neural network methods for transcription unit annotation using precision nuclear run-on (PRO-seq)* <br />
Paul R. Munn, Jay Chia, Charles G. Danko <br />
Unpublished


## Prerequisites

* `Bash >= 4`
* `Python >= 3.5`
* Python modules: `tensorflow >= 2.0, scipy, numpy, matplotlib, getopt`
* `R >= 3.4.2`
* R libraries: `tunits, rqhmm`


## Installation

There is no need for installation of the code in this repository. However, the Tunits R library (written by Andre Martins: [https://github.com/andrelmartins](https://github.com/andrelmartins)) will need to be installed prior to running the HMM code. To install Tunits, follow these steps:
* There are several R packages that will need to be installed to get this working. The first is rqhmm:
* git clone [https://github.com/andrelmartins/QHMM.git](https://github.com/andrelmartins/QHMM.git)
* cd QHMM
* R CMD INSTALL rqhmm
* Then we need to set up the bigwig package, and T-units itself:
* git clone --recursive [https://github.com/andrelmartins/tunits.nhp](https://github.com/andrelmartins/tunits.nhp)
* cd bigWig; R CMD INSTALL bigWig; cd -
* make -C QHMM
* make -C tunits


## Data

The data required to run the programs below can be found in the data directory. This includes the following:
* A pre-trained model produced by the CNN 
* References file for hg19, mm9, and equCab2

Example PRO-seq bigwig files can be found in the GEO, under accession number [GSM1480327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1480327)

## Program usage

### Suggested directory structure

Begin by creating a directory named "data" within the directory you have set up for the python and R scripts. Copy the "models" and "ref_files" directories (and all their sub-directories) into this data directory.

Next, create a "seq" directory under the data directory and copy the bigwig files produced by your PRO-seq assay into this.

Next, create a "bigwigs" directory under the data directory and within this create a "LAB_V2_CNN_V4" directory (this corresponds to the version numbers of the labels and CNN used to build the predictive model). Lastly, under this LAB_V2_CNN_V4 directory, create a directory for each of the chromosomes you will be making predictions for with the following names:

"bigwigs_all_positions_50bp_" + cell type + "_" + chromosome

For example, if you are making predictions for K562 cells, for chromosome 7 you would create a directory with the following name:

"bigwigs_all_positions_50bp_K562_chr7"

So your final directory structure should look like:

```
data
├── bigwigs
│   ├── LAB_V2_CNN_V4
│   │   ├── bigwigs_all_positions_50bp_K562_chr7
│   │   └── bigwigs_all_positions_50bp_K562_chr21
├── models
│   └── LAB_V2_CNN_V4
│   │   └── multiclass-50K-windows-random-center
├── ref_files
│   └── bedbins
└── seq
    └── G1
```

Once this structure is set up, the ```set_up_globals.py``` program need to be edited to point to the data directory you have created. Specifically, line 31 of this code should be changed to:
```
data_folder = '<absolute path above data directory>/data/'
```

For example:
```
data_folder = 'C:/proseq_shapes/data/'
```

### Produce bigwig files of the predictions for each chromosome

```
write-bigwigs-all-positions-50bp.py
```

This code takes as input the chromosome you are making predictions for, the epoch number for the predictive model you are using, the cell type your PRO-seq assay was run on, and the pathways to the bigwigs files generated by your PRO-seq assay (for the plus and minus strands).

Program usage:

```
python write-bigwigs-all-positions-50bp.py \
-c <Chromosome you are making predictions for. Default: chr21> \
-e <Epoch for the model you are using. Default: 3610> \
-l <Cell type for your PRO-seq assay. Default: K562> \
-p <Path for plus bigwig file. Default=seq/G1/G1_plus.bw> \
-m <Path for minus bigwig file. Default=seq/G1/G1_minus.bw>
```

Example:

```
python write-bigwigs-all-positions-50bp.py \
-c chr7 -e 3610 -l K562 -p seq/G1/G1_plus.bw -m seq/G1/G1_minus.bw
```

This program writes the predictions as bigwig files into the appropriate directories under the data/bigwigs directory.

Multiple bigwig files are produced for each 'type' of region that the program is attempting to predict (for example, gene bodies, genes starts, gene ends, etc.), but the Tunits code as written will only use the gene body and gene start pedictions (plus-genebody.bw, minus-genebody.bw, plus-genestart.bw, and minus-genestart.bw). 

Prior to running the Tunits script the gene start files need to be converted to wig files. Additionally, we should filter out predictions of very small read counts (below 0.1 for our purposes) and then merge contiguous regions of read counts so that the Tunits script will run more efficiently.

### Produce thresholded / merged wig files:

Within each of the bigwig directories, run the following comand to convert the bigwig files to wig files:

```
bigWigToWig plus-genestart.bw plus-genestart.wig
bigWigToWig minus-genestart.bw minus-genestart.wig
```

Next, run the following awk script to apply the appropriate threshold:

```
cat plus-genestart.wig | awk 'BEGIN { OFS = "\t" } ; ($4 > 0.1) {print $0}' > plus-genestart-thresh01.wig
cat minus-genestart.wig | awk 'BEGIN { OFS = "\t" } ; ($4 > 0.1) {print $0}' > minus-genestart-thresh01.wig
```

Finally, use bedtools to merge contiguous regions within these thresholded files:

```
bedtools merge -i plus-genestart-thresh01.wig -c 4 -o max > plus-genestart-thresh01-merge.wig
bedtools merge -i minus-genestart-thresh01.wig -c 4 -o max > minus-genestart-thresh01-merge.wig
```

We are now ready to run Tunits' hidden Markov model to smooth the predictions and output the predicted gene bodies as bed files.

### Run Tunits

```
run.tunits.v2.R
```

This code takes as input the chromosome you are making predictions for, the cell type, and the path to your data directory.

Program usage:

```
Rscript --vanilla run.tunits.v2.R <chromosome> <cell type> <path to data directory>
```

Example:

```
Rscript --vanilla run.tunits.v2.R chr7 K562 C:/proseq_shapes/data/
```

Output files:

The program produces one bed file for each strand with the following names:

```
hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma.<cell type>.<chromosome>.3000.preds.plus.bed and
hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma.<cell type>.<chromosome>.3000.preds.minus.bed
```

So, for the example above we would get:

```
hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma.K562.chr7.3000.preds.plus.bed and
hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma.K562.chr7.3000.preds.minus.bed
```



