#!/bin/bash

#AG 13/10/17
#script to run TIS predictions 
#1 Bam to Sam
#2 Matrix of stop to stop regions
#3 Select positive and negative examples
#4 Train models andmake predictions
#5 Filter the predictions
#6 Access prediction accuracy (optional)
#7 Tidy up

usage(){
cat << EOF
usage: $0 options

wrapper to translation initiation site prediction scripts 

OPTIONS:
    -b  bam file of aligned ribo-seq reads
    -g  genome gtf file
    -f  genome fasta file
    -v  validated open reading frames in bed format (optional)
    -o  output directory
    -i  minimum ribo-seq read length (optional)
    -a  maximum ribo-seq read length (optional)
    -h  this help message

example usage: bash prentis.sh -b <riboseq.aligned.bam> -g <genome.gtf> -f <genome.fa> -v <validated_ORFs.bed> -o <output_dir>

EOF
}

while getopts ":b:g:f:i:a:o:v:h" opt; do
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
            echo "-v validated regions in Bed format  $OPTARG"
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

name=$(basename $input_bam)
prefix=${name%.bam}
echo $prefix

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

if [ $minimum_length ] && [ $maximum_length ]; then
    perl scripts/01_make_stop2stop_matrix.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv $minimum_length $maximum_length
else
    perl scripts/01_make_stop2stop_matrix.pl $genome_gtf ${out_dir}/tmp/${prefix}.sam $genome_fasta ${out_dir}/${prefix}_stop2stop.csv
fi

#------------------------------------------------------------------------------------------
#3 Select positive and negative examples
#------------------------------------------------------------------------------------------

if [ $validation_bed ]; then
    perl scripts/02_positive_negative_sets.pl ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/${prefix}_positive.csv ${out_dir}/${prefix}_negative.csv $validation_bed
else
    perl scripts/02_positive_negative_sets.pl ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/${prefix}_positive.csv ${out_dir}/${prefix}_negative.csv
fi

#------------------------------------------------------------------------------------------
#4 Train models andmake predictions
#------------------------------------------------------------------------------------------

if [ ! -d ${out_dir}/model_metrics ]; then
    mkdir ${out_dir}/model_metrics
fi

Rscript scripts/03_run.models.R ${out_dir}/${prefix}_positive.csv ${out_dir}/${prefix}_negative.csv ${out_dir}/${prefix}_stop2stop.csv ${out_dir}/model_metrics/${prefix}_training_glm_summary.csv ${out_dir}/model_metrics/${prefix}_training_glm_coefficients.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_summary.csv ${out_dir}/model_metrics/${prefix}_training_randomforest_variable_importance.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_thresholds_and_metric_scores.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_confusion_matrix.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_summary.csv ${out_dir}/model_metrics/${prefix}_test_randomforest_performance_roc_plot.pdf ${out_dir}/${prefix}_stop2stop_predictions.csv

#------------------------------------------------------------------------------------------
#5 Summarise the predictions
#------------------------------------------------------------------------------------------

perl scripts/04_summarise_predictions.pl $genome_gtf $genome_fasta ${out_dir}/${prefix}_ORF_predictions.bed ${out_dir}/${prefix}_stop2stop_predictions.csv

#------------------------------------------------------------------------------------------
#6 Access prediction accuracy (optional)
#------------------------------------------------------------------------------------------

if [ $validation_bed ]; then
    perl scripts/05_prediction_validation.pl $genome_gtf $genome_fasta $validation_bed ${out_dir}/${prefix}_ORF_predictions.bed > ${out_dir}/${prefix}_report.txt
fi

#------------------------------------------------------------------------------------------
#7 Tidy up
#------------------------------------------------------------------------------------------

rm ${out_dir}/tmp/${prefix}.sam
