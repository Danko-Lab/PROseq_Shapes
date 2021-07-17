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

Much of the data used to run the programs below can be found in the data directory. Raw sequence data and processed files can be found in the GEO, under accession number [GSE130778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130778)

## Program usage

### Figure 2C: Comparison of compartments containing heat shock genes

```
gene_class_compartment_strength_boxplots.py
```

This code takes as input vectors of compartment strengths for compartments containing genes from each of the different groups in our study and produces box plots comparing compartment strength between HS and NHS for genes in these groups.

Program usage:

```
python gene_class_compartment_strength_boxplots.py \
-a <input file for HS HSF-dependent genes> \
-b <input file for NHS HSF-dependent genes> \
-c <input file for HS HSF-independent genes> \
-d <input file for NHS HSF-independent genes> \
-e <input file for HS down regulated genes> \
-f <input file for NHS down regulated genes> \
-g <input file for HS unregulated genes> \
-i <input file for NHS unregulated genes> \
-j <input file for HS silent genes> \
-k <input file for NHS silent genes> \
-o <output file prefix> -l <plot title> -t -h
```

Input parameters:

* Input files: Files containing a single column of numbers representing compartment strengths for each bin along the chromatin for each gene class
* Output file: File name prefix for plot (see output files below)
* -l: Plot title
* -t: Testing flag
* -h: Show usage help

Example:

```
python gene_class_compartment_strength_boxplots.py \
-a data/K562_HS_combined_KRnorm_all_50K_eigen_HSF-dependent_all_genes.txt \
-b data/K562_NHS_combined_KRnorm_all_50K_eigen_HSF-dependent_all_genes.txt \
-c data/K562_HS_combined_KRnorm_all_50K_eigen_HSF-independent_all_genes.txt \
-d data/K562_NHS_combined_KRnorm_all_50K_eigen_HSF-independent_all_genes.txt \
-e data/K562_HS_combined_KRnorm_all_50K_eigen_down_reg.txt \
-f data/K562_NHS_combined_KRnorm_all_50K_eigen_down_reg.txt \
-g data/K562_HS_combined_KRnorm_all_50K_eigen_unreg.txt \
-i data/K562_NHS_combined_KRnorm_all_50K_eigen_unreg.txt \
-j data/K562_HS_combined_KRnorm_all_50K_eigen_silent.txt \
-k data/K562_NHS_combined_KRnorm_all_50K_eigen_silent.txt \
-o output/K562_HS_NHS_gene_class_all_50K_boxplot \
-l 'Comparison of Gene Classes in K562 Cells' \
-t
```

Input text file format:

Single vector of numbers representing compartment strengths for each bin along the chromatin for each gene class.

Output files:

![Comparison of compartments](images/fig_2c.png)



The program produces box plots comparing compartment strength between HS and NHS for compartments containing genes from each group in our study.


### Figure 4B: Virtual 4C tracks

```
virtual_4C_interactions.py 
```

This code takes a file of anchor points and a paired contacts file and produces virtual 4C data in the form of a bedgraph file 

Program usage:

```
python virtual_4C_interactions.py -s <input TSS file> \
-r <input paired contacts file> \
-w <window size> \
-b <bin size> \
-o <output bedgraph file prefix> \
-q <read quality filter> \
-t -h
```

Input parameters:

* Input TSS file: a .bed file containing each of the anchor points you want to examine (the code produces one bedgraph file for each anchor point) - see below for format
* Input paired contacts file: see contact map format below for details
* Window size: Size of region (in bp) around each anchor point
* Bin size: Size of each sub-division of the window (in bp) for which contact end points are summed
* Output bedgraph file: File name prefix for bedgraph file (see output files below)
* -q: Read quality filter (set to zero to include all reads)
* -t: Testing flag
* -h: Show usage help

Example:

```
python virtual_4C_interactions.py \
-s data/HSF1dependently_induced_genes_K562_TSS_only_greater_than_10K.bed \
-r data/K562_HS_combined_merged_nodups.txt \
-o output/K562_HS_combined_virtual_4C_output_10K \
-w 2000000 -b 10000 -q 30 -t
```

Input TSS bed file format:

Standard, 6-column .bed file: Chromosome, Start position of TSS, End position, Name, Score, Strand (note: Start and End should be one bp apart)

Contact map format:

This file consists of paired contacts in the format produced by the [Aiden lab’s Juicer pipeline](https://github.com/aidenlab/juicer) at https://github.com/aidenlab/juicer : 

Column 0: strand 1 (0 for forward, anything else for reverse) <br />
Column 1: chromosome 1 <br />
Column 2: position 1 <br />
Column 3: restriction site fragment 1 (not used by my code) <br />
Column 4: strand 2 (0 for forward, anything else for reverse) <br />
Column 5: chromosome 2 <br />
Column 6: position 2 <br />
Column 7: restriction site fragment 2 (not used by my code) <br />

Note: The merged_nodups.txt file used in the example is too large to store here (approx. 34 Gb) and so will need to be produced by concatenating the fastq files from our Hi-C assays after heat shock (<links>) and running them through the Aiden lab’s Juicer pipeline (see link above).

Output files:

![Virtual 4C plots](images/fig_4b.png)


The program produces one bedgraph file for each anchor point in the input TSS file. These files contain columns for chromosome, start and end positions for each bin, and a read count for each bin.


