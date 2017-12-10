#!/usr/bin/perl -w
use strict;

#to do 05/12/2017
#script to produce ORF bed tracks from randomForest predictions reporting the predicted TIS with the highest score for the longest transcript of each gene

#input files:
my $gtf=$ARGV[0];
my $fasta=$ARGV[1];
my $out_bed=$ARGV[2];
my $matrix_file=$ARGV[3];

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open gtf and get transcript lengths
my %transcripts; #key = gene_id, transcript_id, #value = sum_exon_lengths;

open(GENES1,$gtf) || die "can't open $gtf";      #gft is 1 based
while (<GENES1>){
    unless(/^#/){
        my @b=split("\t");
        my $chr=$b[0];
        my $class=$b[2];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];
        my ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;
        my ($transcript_id) = $b[8] =~ /transcript_id\s"([^\"]+)";/;

        if ($gene_id && $transcript_id){

            if ($class eq "exon"){
                if ($dir eq "+"){
                    for ($start .. $end){
                        $transcripts{$gene_id}{$transcript_id}++;
                    }
                }else{
                    for ($start .. $end){
                        $transcripts{$gene_id}{$transcript_id}++;
                    }
                }
            }
        }
    }
}
close (GENES1);

#select longest transcripts per gene
my %longest_transcript; #key=gene_id, value=transcript_id

for my $gene (keys %transcripts){
    my $longest=0;
    for my $transcript (keys %{ $transcripts{$gene}} ){
        if ($transcripts{$gene}{$transcript} > $longest) {
            $longest_transcript{$gene}=$transcript;
            $longest=$transcripts{$gene}{$transcript};
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#second pass through the genome, find annotated start codons and setup transcript models for longest transcript of each gene

my %gene_start_codon_fwd;
my %gene_stop_codon_fwd;
my %gene_exons_fwd;

my %gene_start_codon_rev;
my %gene_stop_codon_rev;
my %gene_exons_rev;

my %gene_2_chr; #key = gene_id; value = chr

open(GENES2,$gtf) || die "can't open $gtf";      #gft is 1 based
while (<GENES2>){
    unless(/^#/){
        my @b=split("\t");
        my $chr=$b[0];
        my $class=$b[2];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];
        my ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;
        my ($transcript_id) = $b[8] =~ /transcript_id\s"([^\"]+)";/;

        if ($gene_id && $transcript_id){

            #if the transcript is in the list of longest transcripts
            if (exists ( $longest_transcript{$gene_id} )){

                if ($transcript_id eq $longest_transcript{$gene_id}){

                    $gene_2_chr{$gene_id}=$chr;

                    if ($dir eq "+"){ #fwd cases. Use start positions as 5'

                        if ($class eq "start_codon"){
                            $gene_start_codon_fwd{$gene_id}=$start;
                        }
                        if ($class eq "stop_codon"){
                            $gene_stop_codon_fwd{$gene_id}=$start;
                        }
                        if ($class eq "exon"){
                            $gene_exons_fwd{$gene_id}{$start}=$end;
                        }

                    }else{ #revese cases use end as 5'

                        if ($class eq "start_codon"){
                            $gene_start_codon_rev{$gene_id}=$end;
                        }
                        if ($class eq "stop_codon"){
                            $gene_stop_codon_rev{$gene_id}=$end;
                        }
                        if ($class eq "exon"){
                            $gene_exons_rev{$gene_id}{$start}=$end;
                        }
                    }
                }
            }
        }
    }
}
close(GENES2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
##open fasta for codon sequeces
my %fasta_sequences; #key = sequence_name, value=sequence
my $name;
open (FA, $fasta) || die "can't open $fasta";
while (<FA>){
    chomp;
    if (/^>([^\s]+)/){ #take header up to the first space
        $name=$1;
        if ($name =~ /^chr(.*)/){
           $name=$1; #if the chr name have a chr* prefix, remove it  
        }
    }else{
        $fasta_sequences{$name}.=$_;
    }
}
close(FA);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign search regions for each genes
#If there is an annotated leader, we will extend the search region though the leader untill an inframe stop codon is found.

my %stop2stop_stop_coord_fwd; #find the index of the stop codon per gene
my %stop2stop_start_coord_fwd; #find the index of the start codon per gene
my %stop2stop_stop_coord_rev;
my %stop2stop_start_coord_rev;
my %gene_model_fwd;

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;

                if ($_ == $gene_stop_codon_fwd{$gene}){
                    $stop2stop_stop_coord_fwd{$gene}=$model_pos;    #find the index of the stop codon per gene
                }

                if ($_ == $gene_start_codon_fwd{$gene}){
                    $stop2stop_start_coord_fwd{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
            }
        }
    }
}


my %gene_model_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;

                if ($exon_start == $gene_stop_codon_rev{$gene}){
                    $stop2stop_stop_coord_rev{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($exon_start == $gene_start_codon_rev{$gene}){
                    $stop2stop_start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
                $exon_start--;
            }
        }
    }
}

my %start_of_stop2stop_fwd; #parse once to find how far upstream we should go. A.K.A. the beginning of the stop2stop region
for my $gene (keys %stop2stop_start_coord_fwd){

    $start_of_stop2stop_fwd{$gene}=$stop2stop_start_coord_fwd{$gene}; #take the annotateded start codon as default.
    my $chr=$gene_2_chr{$gene};
    my $coord=$stop2stop_start_coord_fwd{$gene}-3;
    my $search=1;

    while ($search){

        if (exists ($gene_model_fwd{$gene}{$coord})){ #search upstream while we are within an annotated leader region.
 
            my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord-0)} -1), 1);     #remember the fasta sequence is zero based
            my $codon2=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord+1)} -1), 1);
            my $codon3=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord+2)} -1), 1);
            my $seq=$codon1.$codon2.$codon3;

            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/){
                $search=0;
            }else{
                $start_of_stop2stop_fwd{$gene}=$coord;
            }
            $coord=$coord-3; #go the next upstream codon.
        }else{
            $search=0;
        }
    }
}

my %start_of_stop2stop_rev;
for my $gene (keys %stop2stop_start_coord_rev){

    $start_of_stop2stop_rev{$gene}=$stop2stop_start_coord_rev{$gene}; #take the annotaed start codon as default.
    my $chr=$gene_2_chr{$gene};
    my $coord=$stop2stop_start_coord_rev{$gene}-3;
    my $search=1;

    while ($search){

        if (exists ($gene_model_rev{$gene}{$coord})){ #search upstream while we are within an annotated leader region.

            my $codon1=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($coord-0) } -1) ,1); #remember the fasta sequence is zero based
            my $codon2=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($coord+1) } -1) ,1);
            my $codon3=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($coord+2) } -1) ,1);
            my $seq=$codon1.$codon2.$codon3;
            $seq=~tr/ACGTacgt/TGCAtgca/;

            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/){
                $search=0;
            }else{
                $start_of_stop2stop_rev{$gene}=$coord;
            }
            $coord=$coord-3; #go the next upstream codon.
        }else{
            $search=0;
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my %gene_score;  #gene = score
my %gene_keeper; #gene = matrix row
#get highest position for each gene

open (MAT, $matrix_file) || die "can't open $matrix_file";      #bed is O based at start, 1 based at end
while (<MAT>){

    my $row=$_;
    chomp($row);

    unless (/^\"\"/){ #skip header
        my @l=split(",",$row);
        my $codon=$l[1];
        my $dir=$l[2];
        my $reads_at_pos=$l[8];
        my $reads_down=$l[9];   
        my $reads_up=$l[10];     
        my $prediction=$l[-2]; #second to last column
        my $pos_prob=$l[-1]; #last column
        $pos_prob =~ tr/"//d;

        my ($CHR, $POS, $DIR, $GENE)= $l[0] =~ /^\"(.*)_(\d+)_(fwd|rev)_([^\"]+)\"$/;

        if ($prediction eq "\"pos\""){
 
            if (exists $gene_score{$GENE}){
                if ($pos_prob > $gene_score{$GENE}){ #new winner         
                    $gene_score{$GENE}=$pos_prob;
                    $gene_keeper{$GENE}=$_;
                }
            }else{ #initalise
                $gene_score{$GENE}=$pos_prob;
                $gene_keeper{$GENE}=$_;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output predicted regions in bed format, classify as tuncated, elongated or matching compared to anottated genes
my $ano=0;
my $tru=0;
my $elo=0;

open (OUT, ">$out_bed") || die "Can't open $out_bed";

#header
print OUT "track_name=annotated_TIS itemRgb=On\n";

for my $gene (keys %gene_keeper){

    my $line=$gene_keeper{$gene};
    chomp($line);
    my @l=split(",",$line);

    my $codon=$l[1];
    my $dir=$l[2];
    my $reads_at_pos=$l[8];
    my $reads_down=$l[9];    
    my $reads_up=$l[10];      
    my $prediction=$l[-2];
    my $pos_prob=$l[-1];
    $pos_prob =~ tr/"//d;

    my ($CHR, $POS, $DIR, $GENE)= $l[0] =~ /^\"(.*)_(\d+)_(fwd|rev)_([^\"]+)\"$/;
 
    #find the nearest in frame stop codon;             
    if ($dir eq "\"fwd\""){ #forward cases
 
        #go though transcript coords, untill we find the stop codon
        my $start_coord=$stop2stop_start_coord_fwd{$GENE};
        my $stop_coord=$stop2stop_stop_coord_fwd{$GENE};
        
        my $model_pos=0;
        my $search=1;
        while ($search){

            if ($model_pos >= $stop_coord){ #if the coord is higher than the stop codon coord
                last;
                $search=0;
            }         

            if ($gene_model_fwd{$GENE}{$model_pos} == $POS){ #the predeicted ORF starts here
                if ($model_pos == $start_coord){
                    print OUT "$CHR\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\tAnnotated\t$pos_prob\t+\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\t0,204,0\n";
                    $ano++;                    
                }
                elsif($model_pos > $start_coord){
                    print OUT "$CHR\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\tTruncation\t$pos_prob\t+\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\t160,160,160\n";
                    $tru++;
                }else{
                    print OUT "$CHR\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\tExtension\t$pos_prob\t+\t".($gene_model_fwd{$GENE}{$model_pos}-1)."\t".($gene_model_fwd{$GENE}{$stop_coord}+2)."\t225,218,0\n";
                    $elo++;      
                }
                $search=0;
            } 
            $model_pos++;
        }

    }else{ #reverse cases

        my $start_coord=$stop2stop_start_coord_rev{$GENE};
        my $stop_coord=$stop2stop_stop_coord_rev{$GENE};

        my $model_pos=0;
        my $search=1;
        while ($search){

            if ($model_pos >= $stop_coord){ #if the coord is higher than the stop codon coord
                $search=0;
            }

            if ($gene_model_rev{$GENE}{$model_pos} == $POS){ #the predicted ORF starts here

                if ($model_pos == $start_coord){
                    print OUT "$CHR\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\tAnnotated\t$pos_prob\t-\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\t0,204,0\n";
                    $ano++;
                }elsif($model_pos > $start_coord){
                    print OUT "$CHR\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\tTruncation\t$pos_prob\t-\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\t160,160,160\n";
                    $tru++;
                }else{
                    print OUT "$CHR\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\tExtension\t$pos_prob\t-\t".($gene_model_rev{$GENE}{$stop_coord}-3)."\t$gene_model_rev{$GENE}{$model_pos}\t225,218,0\n";
                    $elo++;
                }
                $search=0;
            }
            $model_pos++;
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
print "Counts:\n$ano\t#1 Predictions matching annotated start codons\n$tru\t#2 Predicted truncations\n$elo\t#3 Predicted extensions\n";

exit;
