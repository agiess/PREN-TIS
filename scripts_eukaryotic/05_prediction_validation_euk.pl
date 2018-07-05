#!/usr/bin/perl -w
use strict;

#to do 05/12/2017
#script to take a bed file of orf predictions and to report n-termini suppoprt per prediction catagory
#looking for exact matches only

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $nterm=$ARGV[2]; #bed of N-termini predictions
my $pred=$ARGV[3];   #of predictions

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
                            if (exists ($gene_start_codon_fwd{$gene_id})){ #if multiple start codon line take the lower
                                if ($start < $gene_start_codon_fwd{$gene_id}){
                                     $gene_start_codon_fwd{$gene_id}=$start;
                                }
                            }else{
                                $gene_start_codon_fwd{$gene_id}=$start;
                            }
                        }
                        if ($class eq "stop_codon"){
                            if (exists ($gene_stop_codon_fwd{$gene_id})){ #if multiple stop codon line take the lower
                                if ($start < $gene_stop_codon_fwd{$gene_id}){
                                     $gene_stop_codon_fwd{$gene_id}=$start;
                                }
                            }else{
                                $gene_stop_codon_fwd{$gene_id}=$start;
                            }
                        }
                        if ($class eq "exon"){
                            $gene_exons_fwd{$gene_id}{$start}=$end;
                        }

                    }else{ #revese cases use end as 5'

                        if ($class eq "start_codon"){
                            if (exists ($gene_start_codon_rev{$gene_id})){ #if multiple start codon line take the higher
                                if ($end > $gene_start_codon_rev{$gene_id}){
                                     $gene_start_codon_rev{$gene_id}=$end;
                                }
                            }else{
                                $gene_start_codon_rev{$gene_id}=$end;
                            }
                        }
                        if ($class eq "stop_codon"){
                            if (exists ($gene_stop_codon_rev{$gene_id})){ #if multiple stop codon line take the higher
                                if ($end > $gene_stop_codon_rev{$gene_id}){
                                     $gene_stop_codon_rev{$gene_id}=$end;
                                }
                            }else{
                                $gene_stop_codon_rev{$gene_id}=$end;
                            }
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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Open fasta for codon sequeces

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

            my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord-0)} -1), 1); #Fasta sequence is zero based
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

            my $codon1=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($coord-0) } -1) ,1); #Fasta sequence is zero based
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
#open n-termini bed
my %n_term_fwd; #key=chr, key=start_position, value=count
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
            $n_term_fwd{$chr}{$stop-2}=$start; #1st nt in codon
        }else{                                                 #reverse cases
            #five_prime=$stop;
            $n_term_rev{$chr}{$start+2}=$stop; #1st nt in codon
        }
    }
}
close (PEP);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign N-termini to genes
my $total_TIS_count_in_supported_genes=0;
my $supported_ORF_count=0;
my $supported_matching=0;
my $supported_truncation=0;
my $supported_elongation=0;
my %N_terminal_assingment_fwd; #key1 = chr, key2= stop, value = ["Annotated", "Truncated", "Extended"]
my %N_terminal_assingment_rev; #key1 = chr, key2= stop, value = ["Annotated", "Truncated", "Extended"]


for my $gene (keys %gene_model_fwd){
    my $stop_codon_coord = $stop2stop_stop_coord_fwd{$gene};
    my $start_codon_coord = $stop2stop_start_coord_fwd{$gene};  
    my $stop_codon_pos = $gene_model_fwd{$gene}{ $stop_codon_coord };  
	my $chr=$gene_2_chr{$gene};  

    if ($start_codon_coord > 20){ #restrict to genes with at least 21nt of 5' leader, upstream of the start codon (for calculating features in the -20 nt window

        if (exists ($n_term_fwd{$chr}{$stop_codon_pos})){  #if they share a stop codon
            $supported_ORF_count++;
   
            my $start_codon_pos = $gene_model_fwd{$gene}{ $start_codon_coord };
            my $nterm_start = $n_term_fwd{$chr}{$stop_codon_pos};

            if ($nterm_start == $start_codon_pos){
                $supported_matching++;
                $N_terminal_assingment_fwd{$chr}{$stop_codon_pos}="Annotated";
            }elsif( $nterm_start < $start_codon_pos ){
                $supported_truncation++;
                $N_terminal_assingment_fwd{$chr}{$stop_codon_pos}="Truncated";
            }else{
                $supported_elongation++;
                $N_terminal_assingment_fwd{$chr}{$stop_codon_pos}="Extended";
            }
 
            #Count potential TIS's here
            my $search_coord=$start_of_stop2stop_fwd{$gene};

            while ($search_coord < $stop_codon_coord){ #start of stop2stop region to stop codon
 
                my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord-0)} -1), 1);     #remember the fasta sequence is zero based
                my $codon2=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord+1)} -1), 1);
                my $codon3=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord+2)} -1), 1);
                my $seq=$codon1.$codon2.$codon3;

                if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/){
                    $total_TIS_count_in_supported_genes++;
                }
                $search_coord+=3;
            }
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $stop_codon_coord = $stop2stop_stop_coord_rev{$gene};
    my $start_codon_coord = $stop2stop_start_coord_rev{$gene};
    my $stop_codon_pos = $gene_model_rev{$gene}{ $stop_codon_coord };
    my $chr=$gene_2_chr{$gene};

    if ($start_codon_coord > 20){ #restrict to genes with at least 21nt of 5' leader, upstream of the start codon (for calculating features in the -20 nt window

        if (exists ($n_term_rev{$chr}{$stop_codon_pos})){  #if they share a stop codon
            $supported_ORF_count++;

            my $start_codon_pos = $gene_model_rev{$gene}{ $start_codon_coord };
            my $nterm_start = $n_term_rev{$chr}{$stop_codon_pos};

            if ($nterm_start == $start_codon_pos){
                $supported_matching++;
                $N_terminal_assingment_rev{$chr}{$stop_codon_pos}="Annotated";
            }elsif( $nterm_start > $start_codon_pos ){
                $supported_truncation++;
                $N_terminal_assingment_rev{$chr}{$stop_codon_pos}="Truncated";
            }else{
                $supported_elongation++;
                $N_terminal_assingment_rev{$chr}{$stop_codon_pos}="Extended";
            }

            #I can also count poitential TIS's here
            my $search_coord=$start_of_stop2stop_rev{$gene};
            while ($search_coord < $stop_codon_coord){ #start of stop2stop region to stop codon

                my $codon1=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($search_coord-0)} -1), 1);     #remember the fasta sequence is zero based
                my $codon2=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($search_coord+1)} -1), 1);
                my $codon3=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($search_coord+2)} -1), 1);
                my $seq=$codon1.$codon2.$codon3;
                $seq=~tr/ACGTacgt/TGCAtgca/;

                if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/){
                    $total_TIS_count_in_supported_genes++;
                }
                $search_coord+=3;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

            #check if the region shares a stop codon with a peptide
            if (exists ($N_terminal_assingment_fwd{$chr}{$stop-2} )){

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
            
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

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
my $mcc=eval{ ( ($predicted_correct * $correctly_rejected) - ($predicted_incorrect * $supported_not_found) ) / (sqrt ( ($predicted_correct + $predicted_incorrect) * ($predicted_correct + $supported_not_found) * ($correctly_rejected + $predicted_incorrect) * ($correctly_rejected + $supported_not_found) ) ) } || 0;

print "\n";
print "sensitivity\t$sensitivity\n";
print "specificity\t$specificity\n";
print "positive_predictive_value\t$precision\n";
print "matthews_correlation_coefficient\t$mcc\n";

exit;
