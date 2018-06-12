#!/bin/bash

#AG 08/02/18
#script to run TIS predictions 
#1 Bam to Sam
#2 Matrix of stop to stop regions
#3 Select positive and negative examples
#4 Train models and make predictions
#5 Filter the predictions
#6 Access prediction accuracy (optional)
#7 Tidy up

usage(){
cat << EOF
usage: $0 options

wrapper to translation initiation site prediction scripts 

OPTIONS:

    required:
    -b  <file.bam>                 bam file of aligned ribo-seq reads
    -g  <genome.gtf>               genome gtf file
    -f  <genome.fasta>             genome fasta file
    -o  <directory>                output directory

    optional:
    -v  <file.bed>                 validated open reading frames in bed format
    -p  <number between 0 and 1>   the proportion of the 50% most highly expressed genes to use in the positive set (defaults to 1.0 for prokaryotic genomes and 0.1 for eukaryotic genomes)
    -e  No argument                flag for eukaryotic samples (defaults to prokaryotic)
    -d  No argument                flag to enable glm feature selection
    -s  <seed>                     set a seed for repoducibility (defaults to random)
'   -t  <number>                   number of threads for model training/prediction (defaults to 1 thread)
    -i  <number>                   minimum ribo-seq read length (defaults to minimum length of mapped reads in bam file)
    -a  <number>                   maximum ribo-seq read length (defaults to maximum length of mapped reads in bam file)
    -h  No argument                this help message

example usage: bash pren_tis.sh -b <riboseq.aligned.bam> -g <genome.gtf> -f <genome.fa> -v <validated_ORFs.bed> -o <output_dir>

EOF
}

while getopts ":b:g:f:i:a:o:v:p:s:edt:h" opt; do
    case $opt in
        b)
            input_bam=$OPTARG
            echo "-b input bam alignment file $OPTARG"
            ;;
        g)
            genome_gtf=$OPTARG
            echo "-g genome gtf file $OPTARG"
            ;;
        f)
            genome_fasta=$OPTARG
            echo "-f genome fasta file $OPTARG"
            ;;
        i)	
            minimum_length=$OPTARG
            echo "-m minimum ribo-seq read length $OPTARG"
            ;;
        a) 
            maximum_length=$OPTARG
            echo "-i maximum ribo-seq read length $OPTARG"
            ;;
        o) 
            out_dir=$OPTARG
            echo "-o output directory $OPTARG"
            ;;
        v) 
            validation_bed=$OPTARG
            echo "-v validated regions in bed format $OPTARG"
            ;;
        e)
            eukaryotic=1    
            echo "-e eukaryotic flag set"
            ;;
        d)
            use_glm=1
            echo "-d glm feature selection enabled"
            ;;
        p)  
            proportion_of_high_genes=$OPTARG
            ;;
        s)
            seed=$OPTARG
            ;;
        t)  
            threads=$OPTARG
            ;;
        h)
            usage
            exit
            ;;
        ?) 
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done

if [ ! $input_bam ] || [ ! $genome_gtf ] || [ ! $genome_fasta ] || [ ! $out_dir ]; then
    echo "ERROR: something's missing - check your command lines arguments"
    usage
    exit 1
fi

if [ ! -e $input_bam ]; then
    echo "ERROR: could not find the file: $input_bam"
    exit 1
fi

if [ ! -e $genome_gtf ]; then
    echo "ERROR: could not find the file: $genome_gtf"
    exit 1
fi

if [ ! -e $genome_fasta ]; then
    echo "ERROR: could not find the file: $genome_fasta"
    exit 1
fi

if [ ! -d $out_dir ]; then
    echo "ERROR: could not find the directory: $out_dir"
    exit 1
fi

#set default seed
if [ ! $seed ]; then
   seed=$RANDOM;
fi

echo "-s seed set to $seed"

#set default number of threads
if [ ! $threads ]; then
   threads=1
fi 

#set default proportion values
if [ ! $proportion_of_high_genes ]; then

    if [ $eukaryotic ]; then
        proportion_of_high_genes="0.1"
    else
        proportion_of_high_genes="1.0"
    fi

else 

   if (( $(echo "$proportion_of_high_genes <= 0.0" | bc)  )) || (( $(echo "$proportion_of_high_genes > 1.0" | bc) )); then

       echo "ERROR: -p must be between zero and 1"
       exit 1
 
   fi

fi

echo "-t use $threads threads for model training and predictions"
echo "-p use $proportion_of_high_genes of the 50% most highly expressed genes as positive examples"

name=$(basename $input_bam)
prefix=${name%.bam}
echo "Sample: $prefix"

if [ ! -d ${out_dir}/tmp ]; then
    mkdir ${out_dir}/tmp
fi

#------------------------------------------------------------------------------------------
#1 Bam to sam
#------------------------------------------------------------------------------------------

samtools view $input_bam > ${out_dir}/tmp/${prefix}.sam
	
#------------------------------------------------------------------------------------------
#2 Matrix of stop to stop regions
#------------------------------------------------------------------------------------------

if [ $eukaryotic ]; then

    if [ $minimum_length ] && [ $maximum_length ]; then
        perl scripts_eukaryotic/01_make_stop2stop_matrix_euk.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv $minimum_length $maximum_length
    else
        perl scripts_eukaryotic/01_make_stop2stop_matrix_euk.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv
    fi

else

    if [ $minimum_length ] && [ $maximum_length ]; then
        perl scripts/01_make_stop2stop_matrix.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv $minimum_length $maximum_length
    else
        perl scripts/01_make_stop2stop_matrix.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv
    fi

fi

#------------------------------------------------------------------------------------------
#3 Select positive and negative examples
#------------------------------------------------------------------------------------------

if [ $validation_bed ]; then

    perl scripts/02_positive_negative_sets.pl ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/${prefix}_positive_training.csv ${out_dir}/${prefix}_negative_training.csv ${out_dir}/${prefix}_positive_testing.csv ${out_dir}/${prefix}_negative_testing.csv $proportion_of_high_genes $seed $validation_bed

    if [ $? != 0 ] ; then
        echo "An error ocured in selection of the positive and negative sets"
        exit 1;
    fi

else

    perl scripts/02_positive_negative_sets.pl ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/${prefix}_positive_training.csv ${out_dir}/${prefix}_negative_training.csv ${out_dir}/${prefix}_positive_testing.csv ${out_dir}/${prefix}_negative_testing.csv $proportion_of_high_genes $seed

fi

#------------------------------------------------------------------------------------------
#4 Train models and make predictions
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/model_metrics ]; then
    mkdir ${out_dir}/model_metrics
fi
 
if [ $use_glm ]; then

    Rscript scripts/03_run_model_glm.R ${out_dir}/${prefix}_positive_training.csv ${out_dir}/${prefix}_negative_training.csv ${out_dir}/${prefix}_positive_testing.csv ${out_dir}/${prefix}_negative_testing.csv ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/model_metrics/${prefix}_training_glm_summary.csv ${out_dir}/model_metrics/${prefix}_training_glm_coefficients.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_summary.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_thresholds_and_metric_scores.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_confusion_matrix.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_summary.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_roc_plot.pdf ${out_dir}/${prefix}_stop2stop_predictions.csv $threads $seed

    #process GLM coefficients
    perl scripts/03_variable_importance_matrix_GLM.pl ${out_dir}/model_metrics/${prefix}_training_glm_coefficients.csv ${out_dir}/model_metrics/${prefix}_training_glm_coefficients_matrix.csv

    #process RF variable importance
    perl scripts/03_variable_importance_matrix_RF.pl ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance_matrix.csv

    #plot variable importance heatmaps
    Rscript scripts/03_plot_variables_glm.R ${out_dir}/model_metrics/${prefix}_training_glm_coefficients_matrix.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance_matrix.csv ${out_dir}/model_metrics/${prefix}_variable_importance_heatmaps.pdf

else

    Rscript scripts/03_run_model.R ${out_dir}/${prefix}_positive_training.csv ${out_dir}/${prefix}_negative_training.csv ${out_dir}/${prefix}_positive_testing.csv ${out_dir}/${prefix}_negative_testing.csv ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/model_metrics/${prefix}_training_glm_summary.csv ${out_dir}/model_metrics/${prefix}_training_glm_coefficients.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_summary.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_thresholds_and_metric_scores.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_confusion_matrix.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_summary.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_roc_plot.pdf ${out_dir}/${prefix}_stop2stop_predictions.csv $threads $seed

    #process RF variable importance
    perl scripts/03_variable_importance_matrix_RF.pl ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance_matrix.csv

    #plot variable importance heatmaps
    Rscript scripts/03_plot_variables.R ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance_matrix.csv ${out_dir}/model_metrics/${prefix}_variable_importance_heatmaps.pdf

fi

#------------------------------------------------------------------------------------------
#5 Summarise the predictions
#------------------------------------------------------------------------------------------

if [ $eukaryotic ]; then

    perl scripts_eukaryotic/04_summarise_predictions_euk.pl $genome_gtf $genome_fasta ${out_dir}/${prefix}_ORF_predictions.bed ${out_dir}/${prefix}_stop2stop_predictions.csv

else

    perl scripts/04_summarise_predictions.pl $genome_gtf $genome_fasta ${out_dir}/${prefix}_ORF_predictions.bed ${out_dir}/${prefix}_stop2stop_predictions.csv

fi

#------------------------------------------------------------------------------------------
#6 Access prediction accuracy (optional)
#------------------------------------------------------------------------------------------

if [ $validation_bed ]; then

    if [ ! -d ${out_dir}/validation_metrics ]; then
        mkdir ${out_dir}/validation_metrics
    fi

    if [ $eukaryotic ]; then

       perl scripts_eukaryotic/05_prediction_validation_euk.pl $genome_gtf $genome_fasta $validation_bed ${out_dir}/${prefix}_ORF_predictions.bed > ${out_dir}/validation_metrics/${prefix}_performance_on_validated_TIS.txt

    else

       perl scripts/05_prediction_validation.pl $genome_gtf $genome_fasta $validation_bed ${out_dir}/${prefix}_ORF_predictions.bed > ${out_dir}/validation_metrics/${prefix}_performance_on_validated_TIS.txt 

    fi

fi

#------------------------------------------------------------------------------------------
#7 Tidy up
#------------------------------------------------------------------------------------------

rm ${out_dir}/tmp/${prefix}.sam
