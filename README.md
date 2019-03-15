# PREN-TIS

A generalised machine learning framework to exploit experimental ribo-seq read lengths patterns for accurate genome wide identification of translation initiation sites.

## Prerequisites

```
samtools

perl

r

r package: h2o 
r package: dplyr
r package: ggplot2
r package: data.table
```

## Running the program 

```
To launch PREN-TIS simply run the pren_tis.sh wrapper

OPTIONS:
    -b  bam file of aligned ribo-seq reads
    -g  genome gtf file
    -f  genome fasta file
    -o  output directory
    -v  validated open reading frames in bed format (optional)
    -e  flag for eukaryotic samples (defaults to prokaryotic)
    -d  flag to enable glm feature selection (by default no feature selection is used)
    -m  plot heatmaps and read distributions around start codons (defaults to false)
    -p  the proportion of the 50% most highly expressed genes to use in the positive set (defaults to 1.0 for prokaryotic genomes and 0.1 for eukaryotic genomes)
    -s  set a seed for reproducible set selection and model training (defaults to random)
    -t  number of threads for model training/prediction (defaults to 1 thread)
    -i  minimum ribo-seq read length (defaults to minimum length of mapped reads in bam file)
    -a  maximum ribo-seq read length (defaults to maximum length of mapped reads in bam file)
    -h  this help message
```
### Example usage:

```
bash pren_tis.sh -b <aligned.reads.bam> -o <output_directory> -g <genome.gtf> -f <genome.fasta>
```

## Algorithm

The algorithm carries out the following steps:

1) Convert the bam file to a sam file.

2) Calculate values for features at all potential translation sites (TIS), for all genes in the gtf file.

3) Select positive and negative examples for model training.

4) Define training and testing sets.
   (Optionally) Select features with a glm. Alpha and lambda values are tuned with 10 fold cross validation).
   Train a randomforest model using the training set. The number of trees is tuned with 10 fold cross validation. 
   Report the randomforest model metrics on the testing set.
   Score all potential TISs with the randomforest model.

5) Select the most likely predicted TIS for each annotated gene.

6) (optional) Access the predicted TIS against validated TIS from independent methods, such as N-terminal proteomics.

7) (optional) Plot ribo-seq read distirbutions and heatmaps around start codons (defaults to false)
 
8) Remove temporary files.

## Output files

```
prefix_ORF_predictions.bed

model_metrics:
    prefix_training_randomforest_summary.csv 
    prefix_training_randomforest_variable_importance.csv 
    prefix_test_randomforest_performance_thresholds_and_metric_scores.csv 
    prefix_test_randomforest_performance_confusion_matrix.csv 
    prefix_test_randomforest_performance_summary.csv 
    prefix_test_randomforest_performance_roc_plot.pdf 
    prefix_variable_importance_heatmaps.pdf
    prefix_training_glm_summary.csv (optional)
    prefix_training_glm_coefficients.csv (optional)

validation_metrics:
    prefix_performance_on_validated_TIS.txt

heatmaps (optional)
   prefix_heatmaps.pdf 
    
```

## Notes

```
When prokaryotic genomes are supplied the algorithm assumes that exon definitions in the gtf file contain stop codons.

We suggest choosing a value for the proportion of the most highly expressed genes such that 1000-2000 annotated TIS are in the positive set.

This method builds on work published in: [link to paper](https://rdcu.be/brbpQ)
Giess, Adam, Veronique Jonckheere, Elvis Ndah, Katarzyna Chyżyńska, Petra Van Damme, and Eivind Valen. 2017. “Ribosome Signatures Aid Bacterial Translation Initiation Site Identification.” BMC Biology 15 (1): 76.
