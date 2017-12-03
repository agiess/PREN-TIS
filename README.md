# PREN-TIS

A generalised machine learning framework to exploit experimental ribo-seq read lengths patterns for accurate genome wide identification of prokaryotic translation initiation sites.

## Prerequisites

```
samtools

perl

r

r package: h2o 
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
    -e  flag for eukaryotic genomes (defaults to prokaryote)
    -v  validated translation initiation sites in bed format (optional)
    -i  minimum ribo-seq read length (optional)
    -a  maximum ribo-seq read length (optional)
    -h  this help message
```
### Example usage:

```
pren_tis.sh -b <aligned.reads.bam> -o <output_directory> -g <genome.gtf> -f <genome.fasta>
```

## Algorithm

The algorithm carries out the following steps:

1) Convert the bam file to a sam file.

2) Calculate values for features at all potential translation sites (TIS), for all genes in the gtf file.

3) Select positive and negative examples for model training.

4) Define training and testing sets.
   Select features with a glm, using the training set. Alpha and lambda values are tuned with 10 fold cross validation.
   Train a randomforest model using the training set and the selected features. The number of trees is tuned with 10 fold cross validation. 
   Report the randomforest model metrics on the testing set. 
   Score all potential TISs with the randomforest model.

5) Select the most likely predicted TIS for each annotated gene.

6) (optional) Access the predicted TIS against validated TIS from independent methods, such as N-terminal proteomics.

7) Remove temporary files.

## Output files

```
prefix_ORF_predictions.bed

model_metrics:
    prefix_training_glm_summary.csv 
    prefix_training_glm_coefficients.csv 
    prefix_training_randomforest_summary.csv 
    prefix_training_randomforest_variable_importance.csv 
    prefix_test_randomforest_performance_thresholds_and_metric_scores.csv 
    prefix_test_randomforest_performance_confusion_matrix.csv 
    prefix_test_randomforest_performance_summary.csv 
    prefix_test_randomforest_performance_roc_plot.pdf 
    prefix_variable_importance_heatmaps.pdf

validation_metrics:
    prefix_performance_on_validated_TIS.txt
```

## Notes

```
When prokaryotic genomes are supplied the algorithm assumes that exon definitions in the gtf file contain stop codons.
```

