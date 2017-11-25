#!/usr/bin/perl -w
use strict;

#to do 04/10/2017
#script to produce orf bed tracks from randomForest predictions, with probabilities in the bed score
    #find the highest probability for each orf

#input files:
my $gtf=$ARGV[0];
my $fasta=$ARGV[1];
my $out_bed=$ARGV[2];
my $matrix_file=$ARGV[3];

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open gtf file setup genomic hashes, to store start and stop coords per gene, per direction
my %gene_regions_fwd_start;
my %gene_regions_fwd_stop;
my %gene_regions_rev_start;
my %gene_regions_rev_stop;

my %stops_fwd; #chr,pos,gene=1
my %stops_rev; #chr,pos,gene=1;

open(GENES,$gtf) || die "can't open $gtf";        #gtf is 1 based
while (<GENES>){
    unless(/^#/){
        my @b=split("\t");
        my $class=$b[2];
        my $chr=$b[0];
        my $prim5=$b[3];
        my $prim3=$b[4];
        my $dir=$b[6];
        my $gene_id="unknown";
        ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;

        if ($class eq "CDS"){

            if ($dir eq "+"){ #use start positions as begining of feature
                $gene_regions_fwd_start{$gene_id}=$prim5;
                $gene_regions_fwd_stop{$gene_id}=$prim3; #+3 because the cds does not include the stop codon
                $stops_fwd{$chr}{$prim3}=$gene_id;
            }else{
                $gene_regions_rev_start{$gene_id}=$prim3;
                $gene_regions_rev_stop{$gene_id}=$prim5; #-3 becuase the cds does not include the stop codon
                $stops_rev{$chr}{$prim5}=$gene_id;
            }

           # if ($dir eq "+"){ #use start positions as begining of feature
           #     $gene_regions_fwd_start{$gene_id}=$prim5;
           #     $gene_regions_fwd_stop{$gene_id}=$prim3+3; #+3 because the cds does not include the stop codon
           #     $stops_fwd{$chr}{$prim3+3}=$gene_id;
           # }else{
           #     $gene_regions_rev_start{$gene_id}=$prim3;
           #     $gene_regions_rev_stop{$gene_id}=$prim5-3; #-3 becuase the cds does not include the stop codon
           #     $stops_rev{$chr}{$prim5-3}=$gene_id;
           # }

        }
    }
}
close(GENES);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
##open fasta for codon sequeces
my %fasta_sequences; #key = sequence_name, value=sequence
my $name;
open (FA, $fasta) || die "can't open $fasta";
while (<FA>){
    chomp;
    if (/^>([^\s]+)/){ #take header up to the first space
        $name=$1;
    }else{
        $fasta_sequences{$name}.=$_;
    }
}
close(FA);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#group predictrions by stop codon
my %keepers; #chr,dir,pos=predictionRow
my %stopCodonScores; #chr,dir,pos=score

my $rowCount=0;

open (MAT, $matrix_file) || die "can't open $matrix_file";      #bed is O based at start, 1 based at end
while (<MAT>){

    my $row=$_;
    $rowCount++;
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

        my ($CHR, $POS) = $l[0] =~ /^"(.*)_(\d+)_(fwd|rev)"$/;

        if ($dir eq "\"fwd\""){ #forward cases
 
            #search for the next in frame stop codon
            my $search=1;
            my $pos=$POS-1;  #minus one to match zero based fasta index
            while ($search){
                $pos=$pos+3; #only in frame positions
                if ($pos>length($fasta_sequences{$CHR})){ last; } #check that we don't go out of chr limits

                my $seq=substr($fasta_sequences{$CHR},$pos,3);
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                    $search=0;

                    if (exists ($stopCodonScores{$CHR}{$dir}{$pos}) ){
                        if ($pos_prob > $stopCodonScores{$CHR}{$dir}{$pos}){ #new winner
                            $stopCodonScores{$CHR}{$dir}{$pos}=$pos_prob;
                            $keepers{$CHR}{$dir}{$pos}=$row; 
                        }
                    }else{ #initalise
                        $stopCodonScores{$CHR}{$dir}{$pos}=$pos_prob;
                        $keepers{$CHR}{$dir}{$pos}=$row;
                    }
                }
            }  
        }else{ #rev cases
            my $search=1;
            my $pos=$POS-1;
            while ($search){
                $pos=$pos-3;
                if ($pos<=1){ last; } #check that we don't go out of chr limits

                my $seq=reverse(substr($fasta_sequences{$CHR},$pos-2,3));
                $seq=~tr/ACGTacgt/TGCAtgca/;
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                    $search=0;

                    if (exists ($stopCodonScores{$CHR}{$dir}{$pos}) ){
                        if ($pos_prob > $stopCodonScores{$CHR}{$dir}{$pos}){ #new winner
                            $stopCodonScores{$CHR}{$dir}{$pos}=$pos_prob;
                            $keepers{$CHR}{$dir}{$pos}=$row;
                        }
                    }else{ #initalise
                        $stopCodonScores{$CHR}{$dir}{$pos}=$pos_prob;
                        $keepers{$CHR}{$dir}{$pos}=$row;
                    }
                }
            }
        }
    }
}

my $distinctCount=0;
for my $chr (keys %keepers){
    for my $dir (keys %{ $keepers{$chr} } ){
        for my $pos (keys %{ $keepers{$chr}{$dir} } ){
            $distinctCount++;
        }
    }
}
print "there are $rowCount predictions, and $distinctCount distinct predictions\n"; 

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output predicted regions in bed format, classify as tuncated, elongated or matching compared to anottated genes
my $ano=0;
my $tru=0;
my $elo=0;
my $pred_count=0;

open (OUT, ">$out_bed") || die "Can't open $out_bed";

#header
print OUT "track_name=annotated_TIS itemRgb=On\n";

for my $chrH (keys %keepers){

    for my $dirH (keys %{ $keepers{$chrH} } ){

        for my $posH (keys %{ $keepers{$chrH}{$dirH} } ){

            my $line=$keepers{$chrH}{$dirH}{$posH};

            chomp($line);

            unless ($line =~ /^\"\"/){ #skip header

                my @l=split(",",$line);

                my $codon=$l[1];
                my $dir=$l[2];
                my $reads_at_pos=$l[8];
                my $reads_down=$l[9];    
                my $reads_up=$l[10];      
                my $up_down_ratio=$l[11];
                my $down_ratio=$l[12];
                my $up_ratio=$l[13];

                my $prediction=$l[-2];
                my $pos_prob=$l[-1];
                $pos_prob =~ tr/"//d;

                my ($CHR, $POS) = $l[0] =~ /^"(.*)_(\d+)_(fwd|rev)"$/;

                if ($prediction eq "\"pos\""){    
  
                    $pred_count++;
 
                    #find the nearest in frame stop codon;             
                    if ($dir eq "\"fwd\""){ #forward cases
  
                        #get first nucleotide of start codon then proceed upstream 3nt at a time until we find a stop codon (pattern match) 

                        my $search=1;
                        my $limit=length($fasta_sequences{$CHR});  #set upper limit to the end fo the chromosome
                        my $pos=$POS-1;
                        while ($search){
                            $pos=$pos+3;

                            #check that the substing is not smaller than the chr!
                            if ($pos>$limit){ last; }

                            my $seq=substr($fasta_sequences{$CHR},$pos,3);
  
                            #check for stop codon
                            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                                $search=0;

                                #check if the stop codon matches an annotated gene
                                my $pred_start=$POS;
                                my $pred_end=$pos+3; #this should be the end of the stop codon =+2 (+1 for fasta index)

                                if (exists ($stops_fwd{$CHR}{$pred_end} )){

                                    if (exists ($gene_regions_fwd_start{$stops_fwd{$CHR}{$pred_end}})){

                                        my $gene_start=$gene_regions_fwd_start{$stops_fwd{$CHR}{$pred_end}};
  
                                        if ($pred_start==$gene_start){       #exact match
                                            print OUT "$CHR\t".($pred_start-1)."\t$pred_end\tAnnotated\t$pos_prob\t+\t".($pred_start-1)."\t$pred_end\t0,204,0\n";
                                            $ano++;
                                        }elsif($pred_start>$gene_start){     #truncation
                                            print OUT "$CHR\t".($pred_start-1)."\t$pred_end\tTruncation\t$pos_prob\t+\t".($pred_start-1)."\t$pred_end\t160,160,160\n";
                                            $tru++;
                                        }else{                               #extension  
                                            print OUT "$CHR\t".($pred_start-1)."\t$pred_end\tExtension\t$pos_prob\t+\t".($pred_start-1)."\t$pred_end\t225,218,0\n";
                                            $elo++;
                                        } 
                                    }                          
                                }
                            }
                        }  
                    }else{ #reverse cases

                        #get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match)     
                        my $search=1;
                        my $limit=1;  #set upper limit to start of chromosome
                        my $pos=$POS-1;
                        while ($search){
                            $pos=$pos-3;

                            #check that the substing is not bigger or smaller than the chr!
                            if ($pos<=$limit){ last; }

                            my $seq=reverse(substr($fasta_sequences{$CHR},$pos-2,3));
                            $seq=~tr/ACGTacgt/TGCAtgca/;
                            #check for start codon
                            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                                $search=0;

                                #check if the stop codon matches an annotated gene
                                my $pred_start=$POS;  #right most of prediction
                                my $pred_end=$pos-1;  #leftmost of prediction
                         
                                if (exists ($stops_rev{$CHR}{$pred_end} )){

                                    if (exists ($gene_regions_rev_start{$stops_rev{$CHR}{$pred_end}})){

                                        my $gene_start=$gene_regions_rev_start{$stops_rev{$CHR}{$pred_end}};

                                        if ($pred_start == $gene_start){ #exact match
                                            print OUT "$CHR\t".($pred_end-1)."\t$pred_start\tAnnotated\t$pos_prob\t-\t".($pred_end-1)."\t$pred_start\t0,204,0\n";
                                            $ano++;
                                        }elsif($pred_start < $gene_start){ #truncation
                                            print OUT "$CHR\t".($pred_end-1)."\t$pred_start\tTruncation\t$pos_prob\t-\t".($pred_end-1)."\t$pred_start\t160,160,160\n";
                                            $tru++;
                                        }else{ #extension
                                            print OUT "$CHR\t".($pred_end-1)."\t$pred_start\tExtension\t$pos_prob\t-\t".($pred_end-1)."\t$pred_start\t225,218,0\n";
                                            $elo++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#print "Counts:\n$ano\t#1 Predictions matching annotated start codons\n$tru\t#2 Predicted truncations\n$elo\t#3 Predicted extensions\n";

exit;
