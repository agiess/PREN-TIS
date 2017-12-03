#!/usr/bin/perl -w
use strict;

#to do 17/11/2017
#script to produce the matrix for pca, start codon classifiaction
#output is split into lengths.
#for each gene look upstream of the start codon to find the next in frame stop codon

#with fpkm 
#whole genome

#input files:
#gtf exons, start, stop codons
#sam = aligned 
my $gtf=$ARGV[0];
my $sam=$ARGV[1];
my $fasta=$ARGV[2];
my $out_file=$ARGV[3];
my $min_length=$ARGV[4];
my $max_length=$ARGV[5];

my $WINDOW_START=20; #20 nt upstream
my $WINDOW_END=20;   #20 nt downstream

#1) find longest transcripts per gene
#
#2) setup transcript coords (1 .. n) to genomic coords, per gene. accounting for introns
#
#3) setup concaternated trnscript sequences. And a map of transcript sequences to position

#I need to come up with a way to pass disjoint positions to the stop2stop methods
#concatenate strings and store per longest transcript
#store longest transcript positions to genome coords for each gene

#later book keeping
#calculate FPKM and upstream TIS from transcript models + concaternated seq

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
                if ($dir eq "+"){          #I porbably don't need to sepperate fwd+rev here
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
            if (exists ( $transcripts{$gene_id}{$transcript_id} )){

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
close(GENES2);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

print "fasta parsed\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open sam file, store riboseq counts per position
my %counts_fwd; #key = chr, key2 = position, key3 = lenghts, value = counts.
my %counts_rev; #key = chr, key2 = position, key3 = lenghts, value = counts.
my $total_counts=0; #for FPKM 

#If minimum and maxmium read lengths were not provided, calculate them from the sam file 
my $LONGEST_FRAGMENT=1;             #Default_values
my $SHORTEST_FRAGMENT=1000000;      #Default_values

my $find_short=0;
my $find_long=0;

unless ($min_length){
   $find_short=1;
}else{
   $SHORTEST_FRAGMENT=$min_length;
}

unless ($max_length){
   $find_long=1;
}else{
   $LONGEST_FRAGMENT=$max_length;
}

open (SAM, $sam) || die "can't open $sam";      #sam is 1 based
while (<SAM>){
    #skip headers
    unless(/^#/){
        my @b=split("\t");
        my $leftMost=$b[3]; #leftmost position of match 5' for fwd, 3' for rev
        my $flag=$b[1];
        my $chr=$b[2];
        my $mapq=$b[4];
        my $cig=$b[5];
        my $seq=$b[9];
        my $fivePrime;
        my $length=length($seq);

        if ($chr =~ /^chr(.*)/){
           $chr=$1; #if the chr name have a chr* prefix, remove it 
        }

        if ($mapq >= 10){ #mapping uniqness filter
            unless ($flag & 0x4){   #if aligned

                $total_counts++;

                if ($flag & 0x10){  #if rev calculate 5' from the rightmost position

                    #parse cigar for indels and adjust the length of the alignment             
                    my $length5=length($seq);
                    while ($cig =~/(\d+)I/g){   #add to length for insertions
                        $length5+=$1;
                    }
                    while ($cig =~/(\d+)D/g){   #substact from length for deletions
                        $length5-=$1;
                    }

                    $fivePrime=$leftMost+$length5-1;   #sam is 1 based
                    $counts_rev{$chr}{$fivePrime}{$length}++;

                }else{ #if fwd this is easy
                    $fivePrime=$leftMost;
                    $counts_fwd{$chr}{$fivePrime}{$length}++;
                }

                #find the fragment length range for the output matrix 
                if ($find_short){
                    if ($length<$SHORTEST_FRAGMENT){
                        $SHORTEST_FRAGMENT=$length;
                    }
                }

                if ($find_long){
                    if ($length>$LONGEST_FRAGMENT){
                        $LONGEST_FRAGMENT=$length;
                    }
                }
            }
        }else{
            unless ($flag & 0x4){   #if aligned
                $total_counts++;
            }
        }
    }
}
close (SAM);

print "sam parsed\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign search regions for each genes
#If there is an annotoated leader, we will extend the search region though the leader untill an inframe stop codon is found.

my %stop2stop_stop_coord_fwd; #find the index of the stop codon per gene
my %stop2stop_start_coord_fwd; #find the index of the start codon per gene
my %stop2stop_stop_coord_rev; #find the index of the stop codon per gene
my %stop2stop_start_coord_rev; #find the index of the start codon per gene

my %gene_model_fwd;

#gene_model{$gene}{1}=12345
#gene_model{$gene}{2}=12346
#gene_model{$gene}{3}=12347
#gene_model{$gene}{4}=12348
#gene_model{$gene}{5}=12349
#...
#to end of exons

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon
 
        my $model_pos=0;

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)
            
#            print "$gene,$longest_transcript{$gene},$gene_2_chr{$gene},$exon_start,$exon_end\n";

            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;

                if ($_ == $gene_stop_codon_fwd{$gene}){ 
                    $stop2stop_stop_coord_fwd{$gene}=$model_pos;    #find the index of the stop codon per gene
                }

                if ($_ == $gene_start_codon_fwd{$gene}){                    
                    $stop2stop_start_coord_fwd{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;

    #           print "$gene,$longest_transcript{$gene},$gene_2_chr{$gene},$_,$model_pos\n";

            }       
        }
    }
}

my %gene_model_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=($gene_exons_rev{$gene}{$exon_end})+1;

   #why do I need +1 here
   #but not in the fwd example. The gtf shoulød be 1 bases at fwd + rev?

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

#            print "$gene,$longest_transcript{$gene},$gene_2_chr{$gene},$exon_start,$exon_end\n";
           
            while ($exon_start > $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;

                if ($exon_start == $gene_stop_codon_rev{$gene}){
                    $stop2stop_stop_coord_rev{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($exon_start == $gene_start_codon_rev{$gene}){
                    $stop2stop_start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
                $exon_start--;

#                print "$gene,$longest_transcript{$gene},$gene_2_chr{$gene},$exon_start,$model_pos\n";
            }
        }
    }
}            

my %start_of_stop2stop_fwd; #parse once to find how far upstream we should go. A.K.A. the beginning of the stopo2stop region
for my $gene (keys %stop2stop_start_coord_fwd){

    $start_of_stop2stop_fwd{$gene}=$stop2stop_start_coord_fwd{$gene}; #take the annotaed start codon as default.
    my $chr=$gene_2_chr{$gene};
    my $coord=$stop2stop_start_coord_fwd{$gene}-3;    #in some cases thre is no upstream
    my $search=1;    

    while ($search){

        if (exists ($gene_model_fwd{$gene}{$coord})){ #search upstream while we are within an annotated leader region.

            my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord-0)} -1), 1);     #remember the fasta sequence is zero based
            my $codon2=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord+1)} -1), 1);
            my $codon3=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($coord+2)} -1), 1);
            my $seq=$codon1.$codon2.$codon3;
    
#            print "$gene,$chr,$coord,$gene_model_fwd{$gene}{$coord},$seq\n";

            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/){
                $search=0;
#                print "$gene,$longest_transcript{$gene},$chr,$seq,$stop2stop_start_coord_fwd{$gene}\n";
            }else{
                #$start_of_stop2stop_fwd{$gene}=$gene_model_fwd{$gene}{$coord};
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

#            #print "$gene,$chr,$coord,$gene_model_rev{$gene}{$coord},$seq\n";

            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/){
                $search=0;
#                print "$gene,$longest_transcript{$gene},$chr,$seq,$stop2stop_start_coord_rev{$gene}\n";
            }else{
#                $start_of_stop2stop_rev{$gene}=$gene_model_rev{$gene}{$coord};
                $start_of_stop2stop_rev{$gene}=$coord; 
           }
            $coord=$coord-3; #go the next upstream codon.
        }else{
            $search=0;
        }
    }
}


#fwd if greater than start_of_stop2stop_fwd{$gene} and less than ( $stop2stop_stop_coord_fwd{$gene} -3 )
#rev is less than start_of_stop2stop_rev{$gene} and greater than ( $stop2stop_stop_coord_rev{$gene} +3 )

#parse again to setup the whole stop2stop index +/- 20
#for my $gene (keys %gene_model){
#    my $start_index=$stop2stop_start_coord{$gene}
#    my $stop_index= $stop2stop_stop_coord{$gene}

    #take positions betweem start and stop

    #I can possibly do all of the assigning here

#exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#my %start_codon; #key=gene_id, value=pos
#my %gene_info_fwd; #key=gene_id, key2=start, value=stop
#my %gene_info_rev; #key=gene_id, key2=start, value=stop
#my %gene_lengths; #for FPKM

#my %ann_start_codon_fwd; #key=chr, key2=pos, value=1   #used for flags only
#my %ann_start_codon_rev; #key=chr, key2=pos, value=1   
#my %multi_exon;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#summarise and output matrix
open (OUT, ">$out_file") || die;

#header
print OUT "#id,codon,dir,canonical_candidate_sites_ORF_region,annotated_start_site,near_cognate_codon,reads_at_pos,window_reads_downstream,window_reads_upstream,proportion_of_reads_at_position,proportion_of_reads_upstream,proportion_of_reads_downstream,ORF_fpkm";

for my $hpos (-$WINDOW_START .. $WINDOW_END){
    print OUT ",seq_".$hpos;    
}

for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
	for my $pos(-$WINDOW_START .. $WINDOW_END){
         print OUT ",".$pos."_".$len;
	}
}
 
print OUT "\n";


#loop through longest transcripts of forwards genes;
for my $gene (sort keys %start_of_stop2stop_fwd){

    #make suer that we have annotated start and stop codons for this gene
    if (exists ($stop2stop_stop_coord_fwd{$gene})) {

        my $chr=$gene_2_chr{$gene};
        my $start_coord=$start_of_stop2stop_fwd{$gene};
        my $stop_coord=($stop2stop_stop_coord_fwd{$gene})-3;

    #shouldn't this also be +3 (in coord space)
    #might  also fix the -20

        my $start_codon_coords=$stop2stop_start_coord_fwd{$gene};
        my $start_codon_genomic_pos1=$gene_start_codon_fwd{$gene};
        my $start_codon_genomic_pos2=$gene_model_fwd{$gene}{ $start_codon_coords };                   

        while ($start_coord <= $stop_coord){

            $start_coord+=3; #move 3nt at a time inframe with stop codons

            if (($start_coord > 20) && ($start_coord < ($stop_coord-21)) ){  #stay in window bounds (and in frame)
   
                #if the codon is a near cognate, calcuate feaures
                my $codon1=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+0)} -1), 1); #the fasta sequence is zero based
                my $codon2=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+1)} -1), 1);
                my $codon3=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{ ($start_coord+2)} -1), 1);
                my $codon=$codon1.$codon2.$codon3;

                if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){ 

                    my $annotated_start_site="FALSE";
                    my $near_cognate="TRUE";
                    my $reads_at_pos=0;
                    my $window_reads_downstream=0;
                    my $window_reads_upstream=0;
                    my $win_sum=0;
                    my $id=$chr."_".$gene_model_fwd{$gene}{$start_coord}."_fwd"; 

                    #check for start codon
                    if ($start_coord == $stop2stop_start_coord_fwd{$gene}){
                         $annotated_start_site="TRUE";
                    }

                #check up/downstream read proportions        
                #parse all lenghts for count summaries (summmed over all fragment lengths)
                    for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                        for my $win_pos (($start_coord-$WINDOW_START) .. ($start_coord+$WINDOW_END)){
                            if (exists ($counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len})){
                                $win_sum+=$counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len};
                                if ($win_pos==$start_coord){ $reads_at_pos+=$counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len}; }
                                if ($win_pos>$start_coord){ $window_reads_downstream+=$counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len}; }
                                if ($win_pos<$start_coord){ $window_reads_upstream+=$counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len}; }
                            }
                        }
                    }

                    my $proportion_at_position=eval {$reads_at_pos/$win_sum} || 0;
                    my $proportion_upstream=eval {$window_reads_upstream/$win_sum} || 0;
                    my $proportion_downstream=eval {$window_reads_downstream/$win_sum} || 0;

                    #fkpm + codon rank here
                    my ($ORF_FPKM,$codon_rank)=&stop2stop_fwd($chr,$gene,$start_coord);

                    print OUT "$id,$codon,fwd,$codon_rank,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";

                    my $window_seq; #for loop here for coords
                    #Loop through sequence here
                    #output sequence
                    for ($start_coord-$WINDOW_START .. $start_coord+$WINDOW_END){
                        my $nuc=substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{$_}-1), 1);
                        print OUT ",$nuc";
                    }

                    #LOOP THROUGH LENGTHS HERE
                    for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                        #second loop for output
                        for my $win_pos (($start_coord-$WINDOW_START) .. ($start_coord+$WINDOW_END)){
                            my $signal=0;
                            #riboseq signal
                            #%counts; #key = chr, key2 = position, value = counts.
                            if (exists ($counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len} )){
                                $signal+=$counts_fwd{$chr}{ $gene_model_fwd{$gene}{$win_pos} }{$len};
                            }

                            my $fraction_of_window_signal=eval {$signal/$win_sum} || 0;   #this is the proportion of reads per length
                            print OUT ",$fraction_of_window_signal";
                        }
                    }
                    print OUT "\n";
                }
            }
        }
    }
}

#loop through longest transcripts of reverse genes;
for my $gene (sort keys %start_of_stop2stop_rev){

    #make suer that we have annotated start and stop codons for this gene
    if (exists ($stop2stop_stop_coord_rev{$gene})) {

        my $chr=$gene_2_chr{$gene};
        my $start_coord=$start_of_stop2stop_rev{$gene};
        my $stop_coord=($stop2stop_stop_coord_rev{$gene})+3;

        my $start_codon_coords=$stop2stop_start_coord_rev{$gene};
        my $start_codon_genomic_pos1=$gene_start_codon_rev{$gene};
        my $start_codon_genomic_pos2=$gene_model_rev{$gene}{ $start_codon_coords };
        print OUT "$gene,$chr,$start_coord,$stop_coord,$start_codon_coords,$start_codon_genomic_pos1,$start_codon_genomic_pos2\n";

        while ($start_coord <= $stop_coord){

            $start_coord+=3; #move 3nt at a time inframe with stop codons

            if (($start_coord > 20) && ($start_coord < ($stop_coord-21)) ){  #stay in window bounds (and in frame)

                #if the codon is a near cognate, calcuate feaures
                my $codon1=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+0)} -1), 1); #the fasta sequence is zero based
                my $codon2=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+1)} -1), 1);
                my $codon3=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{ ($start_coord+2)} -1), 1);
                my $codon=$codon1.$codon2.$codon3;
                $codon=~tr/ACGTacgt/TGCAtgca/;

                if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){
               
                    my $annotated_start_site="FALSE";
                    my $near_cognate="TRUE";
                    my $reads_at_pos=0;
                    my $window_reads_downstream=0;
                    my $window_reads_upstream=0;
                    my $win_sum=0;
                    my $id=$chr."_".$gene_model_rev{$gene}{$start_coord}."_rev";

                    #check for start codon
                    if ($start_coord == $stop2stop_start_coord_rev{$gene}){
                         $annotated_start_site="TRUE";
                    }

                   #parse all lenghts for count summaries (summmed over all fragment lengths)
                    for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                        for my $win_pos (($start_coord-$WINDOW_START) .. ($start_coord+$WINDOW_END)){
                            if (exists ($counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len})){
                                $win_sum+=$counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len};
                                if ($win_pos==$start_coord){ $reads_at_pos+=$counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len}; }
                                if ($win_pos>$start_coord){ $window_reads_downstream+=$counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len}; }
                                if ($win_pos<$start_coord){ $window_reads_upstream+=$counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len}; }
                            }
                        }
                    }

                    my $proportion_at_position=eval {$reads_at_pos/$win_sum} || 0;
                    my $proportion_upstream=eval {$window_reads_upstream/$win_sum} || 0;
                    my $proportion_downstream=eval {$window_reads_downstream/$win_sum} || 0;

                    #fkpm + codon rank here
                    my ($ORF_FPKM,$codon_rank)=&stop2stop_rev($chr,$gene,$start_coord);

                    print OUT "$id,$codon,rev,$codon_rank,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";

                    my $window_seq; #for loop here for coords
                    #Loop through sequence here
                    #output sequence
                    for ($start_coord-$WINDOW_START .. $start_coord+$WINDOW_END){
                        my $nuc=substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{$_}-1), 1);
                        $nuc=~tr/ACGTacgt/TGCAtgca/;
                        print OUT ",$nuc";
                    }

                    #LOOP THROUGH LENGTHS HERE
                    for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                        #second loop for output
                        for my $win_pos (($start_coord-$WINDOW_START) .. ($start_coord+$WINDOW_END)){
                            my $signal=0;
                            #riboseq signal
                            #%counts; #key = chr, key2 = position, value = counts.
                            if (exists ($counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len} )){
                                $signal+=$counts_rev{$chr}{ $gene_model_rev{$gene}{$win_pos} }{$len};
                            }

                            my $fraction_of_window_signal=eval {$signal/$win_sum} || 0;   #this is the proportion of reads per length
                            print OUT ",$fraction_of_window_signal";
                        }
                    }
                    print OUT "\n";
                }
            }
        }
    }
}

close(OUT);

print "fwd gene processed\nshortest: $SHORTEST_FRAGMENT\nlongest: $LONGEST_FRAGMENT\n";

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub stop2stop_fwd{
    #for a given position (chr + gene_id + coord)
    #calculate fkpm
    #calculate the number of potential start sites upstream of this one

    my $CHR=$_[0];
    my $GENE=$_[1];
    my $TRANSCRIPT_COORD=$_[2];

    my $stop_coord=($stop2stop_stop_coord_fwd{$GENE})-3;
    my $upstream_TIS=0;

    ###
    # Find upstream TIS
    ###

    my $coord=$TRANSCRIPT_COORD;
    while ($coord >= 3){
        $coord-=3;  
         my $codon1=substr($fasta_sequences{$CHR}, (($gene_model_fwd{$GENE}{ ($coord+0)}) -1), 1);
         my $codon2=substr($fasta_sequences{$CHR}, (($gene_model_fwd{$GENE}{ ($coord+1)}) -1), 1);
         my $codon3=substr($fasta_sequences{$CHR}, (($gene_model_fwd{$GENE}{ ($coord+2)}) -1), 1);
         my $codon=$codon1.$codon2.$codon3;

         if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){
             $upstream_TIS++;
         }
    }

    ###
    #find orf FPKM
    ###

    my $orf_sum=0;
    my $orf_length=0;
    for my $current_coord ($TRANSCRIPT_COORD .. $stop_coord){
        $orf_length++;
        for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
            if (exists ($counts_fwd{$CHR}{ $gene_model_fwd{$GENE}{$current_coord} }{$len} )){
                $orf_sum+=$counts_fwd{$CHR}{ $gene_model_fwd{$GENE}{$current_coord} }{$len};
            }
        }
    }

    print "fwd,$GENE,$CHR,$TRANSCRIPT_COORD,$orf_length,$orf_sum\n";

    my $orf_FPKM=eval { (1000000000*$orf_sum)/($total_counts*$orf_length) } || 0;
    return ($orf_FPKM, $upstream_TIS); 
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub stop2stop_rev{

    my $CHR=$_[0];
    my $GENE=$_[1];
    my $TRANSCRIPT_COORD=$_[2];

    my $stop_coord=($stop2stop_stop_coord_rev{$GENE})+3;  #shouldn't this also be -3? if were in  coord space
    my $upstream_TIS=0;

    ###
    #find upstream stop codon
    ###

    my $coord=$TRANSCRIPT_COORD;
    while ($coord >= 3){
        $coord-=3;

        my $codon1=substr($fasta_sequences{$CHR}, (($gene_model_rev{$GENE}{ ($coord+0)}) -1), 1); #the fasta sequence is
        my $codon2=substr($fasta_sequences{$CHR}, (($gene_model_rev{$GENE}{ ($coord+1)}) -1), 1);
        my $codon3=substr($fasta_sequences{$CHR}, (($gene_model_rev{$GENE}{ ($coord+2)}) -1), 1);
        my $codon=$codon1.$codon2.$codon3;
        $codon=~tr/ACGTacgt/TGCAtgca/;

        if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){
            $upstream_TIS++;
        }
    }

    ###
    #find orf FPKM
    ###
  
    my $orf_length=0;
    my $orf_sum=0;
    for my $current_coord ($TRANSCRIPT_COORD .. $stop_coord){
        $orf_length++;
        for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
            if (exists ($counts_rev{$CHR}{ $gene_model_rev{$GENE}{$current_coord} }{$len} )){
                $orf_sum+=$counts_rev{$CHR}{ $gene_model_rev{$GENE}{$current_coord} }{$len};
            }
        }
    }

    print "rev,$GENE,$CHR,$TRANSCRIPT_COORD,$orf_length,$orf_sum\n";

    my $orf_FPKM=eval { (1000000000*$orf_sum)/($total_counts*$orf_length) } || 0;
    return ($orf_FPKM, $upstream_TIS);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
