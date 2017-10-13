PREN-TIS

A generalised machine learning framework to exploit experimental ribo-seq read lengths patterns for accurate genome wide identification of prokaryotic translation initiation sites.

################
# Prerequisites #
################

samtools
perl
r
r package: h2o 
r package: ggplot2

#######################
# Running the program #
#######################

To launch PREN-TIS simply run the pren_tis.sh wrapper

OPTIONS:
    -b  bam file of aligned ribo-seq reads
    -g  genome gtf file
    -f  genome fasta file
    -v  validated open reading frames in bed format (optional)
    -o  output directory
    -i  minimum ribo-seq read length (optional)
    -a  maximum ribo-seq read length (optional)
    -h  this help message

example usage: trim_and_align.sh -f <in.fasta.dir> -o <out_dir> -s Zv9

#############
# Algorithm #
#############

The algorithm carries out the following steps:

1) Convert the bam file to a sam file

2) Calculate values for features at all potential translation sites (TIS), for all genes in the gtf file

3) Select positive and negative examples for model training

4) Define training and testing sets
   Feature selection with a glm using the training set. Alpha and lambda values are tuned with 10 fold cross validation
   Randomforest model training using the training set and the selected features. Number of trees tuned with 10 fold cross validation 
   Report randomforest model metrics on the testing set 
   Score all potential TISs with the randomforest model

5) Select the most likely TIS for each annotated gene

6) (optional) Access the predicted TIS with validated TIS from independent methods, such as N-terminal proteomics

7) Tidy up: delete the temporary sam file

################
# Output files #
################


