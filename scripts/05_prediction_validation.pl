#!/usr/bin/perl -w
use strict;

#to do 09/01/2018
#script to take a bed file of orf predictions and to report n-termini suppoprt per prediction catagory
#looking for exact matches only

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $nterm=$ARGV[2]; #bed of N-termini predictions
my $pred=$ARGV[3];  #bed of predictions

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open fasta 
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
#Open n-termini bed
#Extend n-terminal peptides to downstream in frame stop codons
#Store the start and stop positions of n-terminally supported ORFs

my %n_term_fwd; #key=chr, key=1st nt of stop codon, value=1st nt of start codon
my %n_term_rev;

open (PEP, $nterm) || die "can't open $nterm";      #bed is O based at start, 1 based at end
while (<PEP>){
    unless (/^track/){  #skip header
        my @b=split("\t");
        my $chr=$b[0];
        my $start=$b[1]+1;     #start is zero bases
        my $stop=$b[2];
        my $count=$b[4];
        my $dir=$b[5];

        #assign to metaplots 
        if ($dir eq "+"){                                      #fwd cases
            #five_prime=$start;
            my $stop_position=&closest_inframe_downstream_stop_codon_fwd($chr,$start);
            $n_term_fwd{$chr}{$stop_position}=$start; 
        }else{                                                 #reverse cases
            #five_prime=$stop;
            my $stop_position=&closest_inframe_downstream_stop_codon_rev($chr,$stop);
            $n_term_rev{$chr}{$stop_position}=$stop;
        }
    }
}
close (PEP);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Open genes 
#Find annotated genes that share stop codons with an N-terminal peptide

my $total_TIS_count_in_supported_genes=0;
my $supported_ORF_count=0;
my $supported_matching=0;
my $supported_truncation=0;
my $supported_elongation=0;
my %N_terminal_assingment_fwd; #key1 = chr, key2= stop, value = ["Annotated", "Truncated", "Extended"]
my %N_terminal_assingment_rev; #key1 = chr, key2= stop, value = ["Annotated", "Truncated", "Extended"]

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

            if ($dir eq "+"){ 
 
                if (exists ($n_term_fwd{$chr}{$prim3-2})){ #compare the stop codons

                    $supported_ORF_count++; 

                    #compare start of gene to start of peptide
                    if ($prim5==$n_term_fwd{$chr}{$prim3-2}){
                        $supported_matching++;
                        $N_terminal_assingment_fwd{$chr}{$prim3-2}="Annotated";                 
                    }elsif($prim5 < $n_term_fwd{$chr}{$prim3-2}){
                        $supported_truncation++;
                        $N_terminal_assingment_fwd{$chr}{$prim3-2}="Truncated";
                    }else{
                        $supported_elongation++;
                        $N_terminal_assingment_fwd{$chr}{$prim3-2}="Extended";
                    }


                    #count all TIS
                    my $TIS_count=&count_potential_TIS_in_stop_2_stop_region_fwd($chr,$prim3-2,$prim5);
                                   #chromosome, 1st nt of stop codon, 1st nt of start codon
                    $total_TIS_count_in_supported_genes+=$TIS_count;

                }

            }else{

                if (exists ($n_term_rev{$chr}{$prim5+2})){ #compare the stop codons

                    $supported_ORF_count++;

                    #compare start of gene to start of peptide

                    if ($prim3==$n_term_rev{$chr}{$prim5+2}){
                        $supported_matching++;
                        $N_terminal_assingment_rev{$chr}{$prim5+2}="Annotated"; 
                    }elsif($prim3 > $n_term_rev{$chr}{$prim5+2}){
                        $supported_truncation++;
                        $N_terminal_assingment_rev{$chr}{$prim5+2}="Truncated";
                    }else{
                        $supported_elongation++;
                        $N_terminal_assingment_rev{$chr}{$prim5+2}="Extended";
                    }

                    my $TIS_count=&count_potential_TIS_in_stop_2_stop_region_rev($chr,$prim5+2,$prim3);
                                   #chromosome, 1st nt of stop codon, 1st nt of start codon
                    $total_TIS_count_in_supported_genes+=$TIS_count;
                }
            }
        }
    }
}
close(GENES);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Open predictions. 
#Count how many predictions there are for each gene. And how many of those match the peptide supported ORFs

my $true_positive=0;
my $false_positive=0;
my $total_predictions_in_supported_ORFs=0;
my %n_terminally_supported_ORF_match;

my $positive_annotated=0;
my $positive_truncated=0;
my $positive_extended=0;

open (PRED, $pred) || die "can't open $pred";      #bed is O based at start, 1 based at end
while (<PRED>){
    unless (/^track/){  #skip header
        my @b=split("\t");
        my $chr=$b[0];
        my $start=$b[1]+1;     #start is zero based
        my $stop=$b[2];
        my $dir=$b[5];

        if ($dir eq "+"){

            #check if the prediction shared a stop codon with a N-termnially supported ORF
            if (exists ( $n_term_fwd{$chr}{$stop-2} )){

                $n_terminally_supported_ORF_match{$chr}{$stop-2}=1;

                $total_predictions_in_supported_ORFs++;

                #check if the start codon matches the start of the N-terminal peptide
                if ($start == $n_term_fwd{$chr}{$stop-2}){
                    $true_positive++;

                    if ($N_terminal_assingment_fwd{$chr}{$stop-2} eq "Annotated"){ $positive_annotated++; 
                    }elsif ($N_terminal_assingment_fwd{$chr}{$stop-2} eq "Truncated"){ $positive_truncated++;   
                    }elsif ($N_terminal_assingment_fwd{$chr}{$stop-2} eq "Extended"){ $positive_extended++; }
    
                }else{
                    $false_positive++;
                }
            }
        }else{ #reverse cases

            #check if the prediction shared a stop codon with a N-termnially supported ORF
            if (exists ( $n_term_rev{$chr}{$start+2} )){

                $n_terminally_supported_ORF_match{$chr}{$start+2}=1;

                $total_predictions_in_supported_ORFs++;

                #check if the start codon matches the start of the N-terminal peptide
                if ($stop == $n_term_rev{$chr}{$start+2}){
                    $true_positive++;

                    if ($N_terminal_assingment_rev{$chr}{$start+2} eq "Annotated"){ $positive_annotated++;
                    }elsif ($N_terminal_assingment_rev{$chr}{$start+2} eq "Truncated"){ $positive_truncated++;                
                    }elsif ($N_terminal_assingment_rev{$chr}{$start+2} eq "Extended"){ $positive_extended++; }

                }else{
                    $false_positive++;
                }
            }
        }
    }
}

close (PRED);

my $ORFs_found_count=0;

for my $chr (keys %n_terminally_supported_ORF_match){
    for (keys %{ $n_terminally_supported_ORF_match{$chr} }) {
        $ORFs_found_count++;
    }
}

#TRUE_POSITIVE  = The number of predicted ORFs that agree with N-terminal peptides
#FALSE_POSITIVE = The number of predicted ORFs that disagree with N-terminal peptide
#TRUE_NEGATIVE  = The number of potential start sites in ORFs with N-terminal support, that were neither predicted nor supported
#FALSE_NEGATIVE = The number of ORFs with N-terminal support that were not predicted

my $false_negative=$supported_ORF_count-$ORFs_found_count;
my $true_negative=$total_TIS_count_in_supported_genes-$supported_ORF_count-$false_positive;

print "TOTAL_PREDICTIONS\t$total_predictions_in_supported_ORFs\n";
print "TOTAL_SUPPORTED_ORFS\t$supported_ORF_count\n";
print "\n";

print "TRUE_POSITIVE\t$true_positive\n";
print "FALSE_POSITIVE\t$false_positive\n";
print "TRUE_NEGATIVE\t$true_negative\n";
print "FALSE_NEGATIVE\t$false_negative\n";

print "\nOf the true positives:\n";
print "Supported_matching:\t$positive_annotated of $supported_matching sites predicted\n";
print "Supported_truncation:\t$positive_truncated of $supported_truncation sites predicted\n";
print "Supported_elongation:\t$positive_extended of $supported_elongation sites predicted\n";

my $predicted_correct=$true_positive;
my $predicted_incorrect=$false_positive;
my $supported_not_found=$false_negative;
my $correctly_rejected=$true_negative;

my $sensitivity=eval{$predicted_correct/($predicted_correct+$supported_not_found)} || 0 ;
my $specificity=eval{$correctly_rejected/($predicted_incorrect+$correctly_rejected)} || 0;
my $precision=eval{$predicted_correct/($predicted_correct+$predicted_incorrect)} || 0;;

print "\n";
print "sensitivity\t$sensitivity\n";
print "specificity\t$specificity\n";
print "positive_predictive_value\t$precision\n";

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub closest_inframe_downstream_stop_codon_fwd{

    my $chr=$_[0];
    my $start_position=$_[1];

    if (exists ($fasta_sequences{$chr})){

        my $search=1;
        while ($search){ 
            $start_position=$start_position+3;
                    
            #check that the position is not larger than the chromosome
            if ($start_position>((length($fasta_sequences{$chr})-1))){ last; }

            my $seq=substr($fasta_sequences{$chr},($start_position-1),3);
             
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
        }
    }    
    return ($start_position); #the 1st nuclotide of the stop codon
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub closest_inframe_downstream_stop_codon_rev{
    
    my $chr=$_[0];
    my $start_position=$_[1];
    
    if (exists ($fasta_sequences{$chr})){
        
        my $search=1; 
        while ($search){  #find upstream inframe stop codon
            $start_position=$start_position-3;
                    
            #check that the substing is not smaller than the chromosome
            if ($start_position<0){ last; }
            
            my $seq=reverse(substr($fasta_sequences{$chr},$start_position-3,3));
            $seq=~tr/ACGTacgt/TGCAtgca/;
                    
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
        }

    }
    return $start_position; #the 1st nuclotide of the stop codon
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub count_potential_TIS_in_stop_2_stop_region_fwd{

    my $chr=$_[0];
    my $stop_codon=$_[1];
    my $start_codon=$_[2];

    my $total_potential_TIS=0;

    if (exists ($fasta_sequences{$chr})){

        #start at stop codon and procees upstream 3nt at a time until we find a stop codon.
        #or have gone more than 999nt past the start codon.
    
        my $search=1;
        my $count=0;
        my $upstream_of_TIS=0;
 
        while ($search){ 
            $stop_codon=$stop_codon-3;

            #check that the substing is not smaller than the chr!
            if ($stop_codon<0){ last; }
 
            my $seq=substr($fasta_sequences{$chr},($stop_codon-1),3);

            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
    
            if ($count>=999){ $search=0; }
          
            if ($upstream_of_TIS) { $count+=3; }

            if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ $total_potential_TIS++; }

            if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }

        }
    }
    return($total_potential_TIS);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub count_potential_TIS_in_stop_2_stop_region_rev{

    my $chr=$_[0];
    my $stop_codon=$_[1];
    my $start_codon=$_[2];

    my $total_potential_TIS=0;

    if (exists ($fasta_sequences{$chr})){

        #start at stop codon and procees upstream 3nt at a time until we find a stop codon.
        #or have gone more than 999nt past the start codon.
    
        my $search=1;
        my $count=0;  
        my $upstream_of_TIS=0;
 
        while ($search){ 
            $stop_codon=$stop_codon+3;
 
            #check that the substing is not larger than the chr
            if ($stop_codon>=(length($fasta_sequences{$chr})-1)){ last; }
 
            my $seq=reverse(substr($fasta_sequences{$chr},$stop_codon-3,3));
            $seq=~tr/ACGTacgt/TGCAtgca/;
 
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
            
            if ($count>=999){ $search=0; }
            
            if ($upstream_of_TIS) { $count+=3; }

            if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ $total_potential_TIS++; }

            if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }

        }
    }
    return($total_potential_TIS);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
