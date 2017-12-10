#!/usr/bin/perl -w
use strict;

#to do 05/12/2017
#script to take a bed file of orf predictions and to report n-termini suppoprt per prediction catagory
#looking for exact matches only

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $nterm=$ARGV[2]; #bed of N-termini predictions
my $bed=$ARGV[3];   #of predictions

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
my %supported_genes; #later we want to calculate on only gene with support
my $n_sum_ano=0;
my $n_sum_tru=0;
my $n_sum_elo=0;
my %n_term_assignment_fwd; #chr, start, stop = ann/tru/elo
my %n_term_assignment_rev; #chr, stop, start = ann/tru/elo
my $total_TIS=0;; #count all of the potential TIS is supported genes

for my $gene (keys %gene_model_fwd){
    my $stop_codon_coord = $stop2stop_stop_coord_fwd{$gene};
    my $start_codon_coord = $stop2stop_start_coord_fwd{$gene};  
    my $stop_codon_pos = $gene_model_fwd{$gene}{ $stop_codon_coord };  
	my $chr=$gene_2_chr{$gene};  

    if ($start_codon_coord > 20){ #restrict to genes with at least 21nt of 5' leader, upstream of the start codon (for calculating features in the -20 nt window

        if (exists ($n_term_fwd{$chr}{$stop_codon_pos})){  #if they share a stop codon
            $supported_genes{$gene}=1;

            my $start_codon_pos = $gene_model_fwd{$gene}{ $start_codon_coord };
            my $nterm_start = $n_term_fwd{$chr}{$stop_codon_pos};

            if ($nterm_start == $start_codon_pos){
                $n_sum_ano++;
                $n_term_assignment_fwd{$chr}{$stop_codon_pos}{$nterm_start}="ann";
            }elsif( $nterm_start < $start_codon_pos ){
                $n_sum_tru++;
                $n_term_assignment_fwd{$chr}{$stop_codon_pos}{$nterm_start}="tru";
            }else{
                $n_sum_elo++; 
                $n_term_assignment_fwd{$chr}{$stop_codon_pos}{$nterm_start}="elo";
            }
 
            #I can also count potential TIS's here
            my $search_coord=$start_of_stop2stop_fwd{$gene};

            while ($search_coord < $stop_codon_coord){ #start of stop2stop region to stop codon
 
                my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord-0)} -1), 1);     #remember the fasta sequence is zero based
                my $codon2=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord+1)} -1), 1);
                my $codon3=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($search_coord+2)} -1), 1);
                my $seq=$codon1.$codon2.$codon3;

                if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/){
                    $total_TIS++;
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
            $supported_genes{$gene}=1;

            my $start_codon_pos = $gene_model_rev{$gene}{ $start_codon_coord };
            my $nterm_start = $n_term_rev{$chr}{$stop_codon_pos};

            if ($nterm_start == $start_codon_pos){
                $n_sum_ano++;
                $n_term_assignment_rev{$chr}{$stop_codon_pos}{$nterm_start}="ann";
            }elsif( $nterm_start > $start_codon_pos ){
                $n_sum_tru++;
                $n_term_assignment_rev{$chr}{$stop_codon_pos}{$nterm_start}="tru";
            }else{
                $n_sum_elo++;
                $n_term_assignment_rev{$chr}{$stop_codon_pos}{$nterm_start}="elo";
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
                    $total_TIS++;
                }
                $search_coord+=3;
            }
        }
    }
}

##got to here
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $count_supported=keys %supported_genes;

#print "matching annotation:  $n_sum_ano\n";
#print "predicted truncation: $n_sum_tru\n";
#print "predoicted extension: $n_sum_elo\n";
#print "supported TIS $count_supported\n";
#print "potential TIS $total_TIS\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open bed of predictions
my %call; 
my $pred_sum_ano=0;
my $pred_sum_tru=0;
my $pred_sum_elo=0;

open (BED, $bed) || die;
while (<BED>){

    unless (/^track/){  #skip header
        my @b=split("\t");
        my $chr=$b[0];
        my $start=$b[1]+1;     #start is zero based
        my $stop=$b[2];
        my $count=$b[4];
        my $dir=$b[5];
        my $type=$b[3]; #Annotated Truncation Extension

        if ($type eq "Annotated"){  $pred_sum_ano++; }
        if ($type eq "Extension"){  $pred_sum_elo++; }
        if ($type eq "Truncation"){ $pred_sum_tru++; }

        if ($dir eq "+"){

            $stop=$stop-2; #1st nt of stop codon

            #check if the region shares a stop codon with a peptide
            if (exists ($n_term_assignment_fwd{$chr}{$stop}) ){  

                for my $nstart (keys %{ $n_term_assignment_fwd{$chr}{$stop} } ){

                    if ($type eq "Annotated"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){
                            if ($start == $nstart){ $call{$chr}{$stop}="ano_ano"; }
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){ 
                            $call{$chr}{$stop}="ano_tru"; 
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            $call{$chr}{$stop}="ano_elo";
                        }

                    #check for extensions
                    }elsif ( $type eq "Extension"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){ 
                            $call{$chr}{$stop}="elo_ano";
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){ 
                            $call{$chr}{$stop}="elo_tru"; 
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            if ($start == $nstart){ $call{$chr}{$stop}="elo_elo"; }
                        } 

                    #otherwise truncation
                    }elsif ( $type eq "Truncation"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){ 
                            $call{$chr}{$stop}="tru_ano";
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){
                            if ($start == $nstart){ $call{$chr}{$stop}="tru_tru"; }
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            $call{$chr}{$stop}="tru_elo";
                        }
                    }
                }
            }
        #rev
        }else{

            #for rev, start is 3' (stop)
            #stop is 5' (start)
            $start=$start+2;

            #check if the region shares a stop codon with a peptide
            if (exists ($n_term_assignment_rev{$chr}{$start} ) ){

                for my $nstop (keys %{ $n_term_assignment_rev{$chr}{$start}} ){

                    if ($type eq "Annotated"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ann"){ 
                            if ($stop == $nstop){ $call{$chr}{$start}="ano_ano"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){ 
                            $call{$chr}{$start}="ano_tru";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            $call{$chr}{$start}="ano_elo";
                        }

                    }elsif ($type eq "Extension"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            if ($stop == $nstop){ $call{$chr}{$start}="elo_elo"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){ 
                            $call{$chr}{$start}="elo_tru";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ano"){ 
                            $call{$chr}{$start}="elo_ano"; 
                        }
                 
                    }elsif ($type eq "Truncation"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){
                            if ($stop == $nstop){ $call{$chr}{$start}="tru_tru"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ann"){ 
                            $call{$chr}{$start}="tru_ano";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            $call{$chr}{$start}="tru_elo"; 
                        }
                    }
                }
            }   
        }
    }
}
close(BED);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $ano_ano=0;
my $ano_tru=0;
my $ano_elo=0;
my $tru_ano=0;
my $tru_tru=0;
my $tru_elo=0;
my $elo_ano=0;
my $elo_tru=0;
my $elo_elo=0;

for my $c (keys %call){
    for my $p (keys %{$call{$c}}){
        if ($call{$c}{$p} eq "ano_ano"){ $ano_ano++;}
        if ($call{$c}{$p} eq "ano_tru"){ $ano_tru++;}
        if ($call{$c}{$p} eq "ano_elo"){ $ano_elo++;}
        if ($call{$c}{$p} eq "tru_ano"){ $tru_ano++;}
        if ($call{$c}{$p} eq "tru_tru"){ $tru_tru++;}
        if ($call{$c}{$p} eq "tru_elo"){ $tru_elo++;}
        if ($call{$c}{$p} eq "elo_ano"){ $elo_ano++;}
        if ($call{$c}{$p} eq "elo_tru"){ $elo_tru++;}
        if ($call{$c}{$p} eq "elo_elo"){ $elo_elo++;}
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
print ",,nterminal_support,,\n";
print "predictions,,ann,elo,tru\n";
print ",sum,$n_sum_ano,$n_sum_elo,$n_sum_tru\n";
print "ann,$pred_sum_ano,$ano_ano,$ano_elo,$ano_tru\n";
print "ext,$pred_sum_elo,$elo_ano,$elo_elo,$elo_tru\n";
print "tru,$pred_sum_tru,$tru_ano,$tru_elo,$tru_tru\n";

#summary stats
#True +ve: Supported and predicted ORFs 
#False +ve: Predicted ORFs that disagreed with supported positions.
#False -ve: Supported genes where no ORF was predicted.
#True -ve:  All potential TIS in supported stop2stop regions, that were neither predicted nor supported.

print "\nthere are $count_supported, genes with n-terminal info\n";
print "there are $total_TIS, candidate TIS in all supported genes\n";

my $predicted_correct=$ano_ano+$elo_elo+$tru_tru;
my $predicted_incorrect=$ano_elo+$ano_tru+$elo_ano+$elo_tru+$tru_ano+$tru_elo;
my $supported_not_found=$count_supported-$predicted_correct-$predicted_incorrect;
my $correctly_rejected=$total_TIS-$predicted_correct-$predicted_incorrect-$supported_not_found;

my $sensitivity=$predicted_correct/($predicted_correct+$supported_not_found);
my $specificity=$correctly_rejected/($predicted_incorrect+$correctly_rejected);
my $precision=$predicted_correct/($predicted_correct+$predicted_incorrect);

print "\ntrue positive: $predicted_correct\n";
print "true negative: $correctly_rejected\n";
print "false positive: $predicted_incorrect\n";
print "false negative: $supported_not_found\n";

print "\nsensitivity: $sensitivity\n";
print "specificity: $specificity\n";
print "positive_predictive_value: $precision\n";

exit;
