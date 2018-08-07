#!/usr/bin/perl -w
use strict;

#06/08/2018
#script to count 5' and 3' fragment lengths around start codons

#input files:
my $gtf=$ARGV[0];
my $bam=$ARGV[1];
my $outDir=$ARGV[2];

my ($prefix)=$bam=~/([^\/]+).bam$/;

my $START_UPSTREAM=-100; 
my $START_DOWNSTREAM=500;  

my $PERCENTAGE_TO_EXCLUDE=0.1;

#exclude genes with leaders shorter than START_UPSTREAM
#exclude genes with CDSs shorter than START_DOWNSTREAM 

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
        my $gene_biotype="NA";
        if ($b[8] =~ /gene_biotype\s"([^\"]+)";/){
            $gene_biotype=$1;
        }

        #strip the chr prefix;
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

        if ($gene_id && $transcript_id){

            if ($gene_biotype eq "protein_coding"){ #restrict to protien coding genes (control for ncRNAs)

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
}
close (GENES1);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

        #strip the chr prefixc
        if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
            $chr=$1;
        }

        if ($gene_id && $transcript_id){

            if (exists ( $longest_transcript{$gene_id} )){    #if the transcript is in the list of longest transcripts

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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup trascript models

my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %five_prime_most_coord_fwd;
my %three_prime_most_coord_fwd;

#gene_model{$gene}{0}=12345   #1st nt of start codon
#gene_model{$gene}{1}=12346
#gene_model{$gene}{2}=12347
#gene_model{$gene}{3}=12348
#gene_model{$gene}{4}=12349
#...
#to end of exons             #last nt of stop codon

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;
        $five_prime_most_coord_fwd{$gene}=$model_pos;  #initalise the 5' to the first coord coord

        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};

            #fwd exons are in ascending order
            # start(-1)-> 100958 100975
            #             101077 101715 <-end(+1)

            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;

                if ($_ == $gene_stop_codon_fwd{$gene}){
                    $stop_coord_fwd{$gene}=$model_pos;    #find the index of the stop codon per gene
                }

                if ($_ == $gene_start_codon_fwd{$gene}){
                    $start_coord_fwd{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
            }
        }
        $three_prime_most_coord_fwd{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

my %gene_model_rev;
my %start_coord_rev;
my %stop_coord_rev;

#5' is the #gene_model{$gene}{0}
#3' is the last coord
my %five_prime_most_coord_rev;
my %three_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon

        my $model_pos=0;
        $five_prime_most_coord_rev{$gene}=$model_pos;  #initalise the 5' to the first coord coord

        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};

            #rev exons are sorted in decending order  
            #           447087 447794 <-start(+1)
            # end(-1)-> 446060 446254

            while ($exon_start >= $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;

                if ($exon_start == $gene_stop_codon_rev{$gene}){
                    $stop_coord_rev{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($exon_start == $gene_start_codon_rev{$gene}){
                    $start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                $model_pos++;
                $exon_start--;
            }
        }
        $three_prime_most_coord_rev{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Transcript coord key
#Transcript 5' coord   0                                  #1st nt of annotated transcript
#Transcript 3' coord   $three_prime_most_coord_???{$gene} #last nt of transcript
#Start codon coord     $start_coord_???{$gene}            #1st nt in start codon
#Stop codon coord      $stop_coord_???{$gene}             #1st nt in stop codon
#$gene_model_fwd{$gene}{$coord}==genomic position

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#loop though genes and assign leader, CDS and trailer regions to hashes for quick searching. (added start + stop codons).
my %start_codons_search_fwd; #key1 = chr, key2 = position, value = gene_id(s) ; window position
my %start_codons_search_rev; 

#for initalising
my %start_codon_meta_positions_5; #key1 = gene_id, key2=metapos, key3=value
my %start_codon_meta_positions_3;

my %start_region_signal_5; #gene_id = sum signal_at_start_region
my %start_region_signal_3; 

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};

    my $include_leader=0;
 
    #check if there is 100nt of leader.
    if (exists ($gene_model_fwd{$gene}{($start_coord+$START_UPSTREAM)})){  #-100
        if (($start_coord+$START_UPSTREAM) >= $five_prime_most_coord_fwd{$gene}) { #restrict to genes with leader of a suitable length
             $include_leader=1;
        }
    }

    #should I also exclude CDS's shorter than 500 nt?
    if ($stop_coord-$start_coord < $START_DOWNSTREAM){
        $include_leader=0;
    }

    if ($include_leader){
        for my $pos ($START_UPSTREAM .. $START_DOWNSTREAM){   #-100 to 500
          
            $start_codon_meta_positions_5{$gene}{$pos}=0;
            $start_codon_meta_positions_3{$gene}{$pos}=0;
            $start_region_signal_5{$gene}=0;
            $start_region_signal_3{$gene}=0;
            my $genomic_coord=$gene_model_fwd{$gene}{($start_coord+$pos)};
            if (exists ($start_codons_search_fwd{$chr}{$genomic_coord})){
                $start_codons_search_fwd{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                $start_codons_search_fwd{$chr}{$genomic_coord}=$gene.";".$pos;
            }
        }
    }
}

for my $gene (keys %gene_model_rev){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};

    my $include_leader=0;

    #check if there is 100nt of leader
    if (exists ($gene_model_rev{$gene}{($start_coord+$START_UPSTREAM)})){  #-100
        if (($start_coord + $START_UPSTREAM) >= $five_prime_most_coord_rev{$gene}){ #restrict to genes with leader of a suitable length
            $include_leader=1;
        }
    }

    #check the CDS is long enough to plot
    if ($stop_coord-$start_coord < $START_DOWNSTREAM){
        $include_leader=0;
    }

    if ($include_leader){
        for my $pos ($START_UPSTREAM .. $START_DOWNSTREAM){   #-100 to 500

            $start_codon_meta_positions_5{$gene}{$pos}=0;
            $start_codon_meta_positions_3{$gene}{$pos}=0;
            $start_region_signal_5{$gene}=0;
            $start_region_signal_3{$gene}=0;

            my $genomic_coord=$gene_model_rev{$gene}{($start_coord+$pos)};

            if (exists ($start_codons_search_rev{$chr}{$genomic_coord})){
                $start_codons_search_rev{$chr}{$genomic_coord}.=",".$gene.";".$pos;  #linker == ";"
            }else{
                $start_codons_search_rev{$chr}{$genomic_coord}=$gene.";".$pos;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open and assign TCPseq counts
my %start_codon_meta_positions_by_length_5; #gene, meta position (-50 to 500), length = count
my %start_codon_meta_positions_by_length_3;

open BAM,"samtools view $bam |";
while(<BAM>){

    next if(/^(\@)/);        #skipping the header lines
    s/\n//;  s/\r//;         #remove new line
    my @sam = split(/\t+/);  #split SAM line into array

    my $leftMost=$sam[3]; #leftmost position of match 5' for fwd, 3' for rev
    my $flag=$sam[1];
    my $chr=$sam[2];
    my $mapq=$sam[4];
    my $cigar=$sam[5];
    my $seq=$sam[9];
    my $threePrime;
    my $fivePrime;

    if ($chr=~/chr(.*)/){ #ensembl chromosome names do not contain the "chr" prefix
        $chr=$1;
    }

    unless ($flag & 0x4){   #if aligned

        if ($mapq >= 10){     #mapping uniqnes filter

            if ($flag & 0x10){  #Reverse reads. Starting from the leftmost position parse the cigar and check if matching positions overlap leaders or cds's

                #assign 5' amd 3' to positions
                $threePrime=$leftMost;

                #parse cigar for indels and adjust the length of the alignment             
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){   #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){   #substact from length for deletions
                    $length-=$1;
                }
                $fivePrime=$leftMost+($length-1);              #SAM is 1 based

                #assign to metaplots
                if (exists ($start_codons_search_rev{$chr}{$threePrime})){
                    my @over1=split(",",$start_codons_search_rev{$chr}{$threePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $start_region_signal_3{$gene}++;
                        $start_codon_meta_positions_3{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

                #same again for the 5 prime
                if (exists ($start_codons_search_rev{$chr}{$fivePrime})){
                    my @over1=split(",",$start_codons_search_rev{$chr}{$fivePrime});
                    for (@over1){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+)$/;  #/(\w+)_(-?\w+)/;
                        $start_region_signal_5{$gene}++;
                        $start_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }

            }else{ #if fwd 3' == sam coordinate (leftmost) + read length 

                #parse cigar for indels and adjust the length of the alignment             
                my $length=length($seq);
                while ($cigar =~/(\d+)I/g){   #add to length for insertions
                    $length+=$1;
                }
                while ($cigar =~/(\d+)D/g){   #substact from length for deletions
                    $length-=$1;
                }

                $threePrime=$leftMost+($length-1);              #SAM is 1 based
                $fivePrime=$leftMost;

                #assign 3' to metaplots
                if (exists ($start_codons_search_fwd{$chr}{$threePrime})){
                    my @over4=split(",",$start_codons_search_fwd{$chr}{$threePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $start_region_signal_3{$gene}++;
                        $start_codon_meta_positions_3{$gene}{$meta_pos}+=1;    
                        $start_codon_meta_positions_by_length_3{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }
 
                #do that again for five prime
                if (exists ($start_codons_search_fwd{$chr}{$fivePrime})){
                    my @over4=split(",",$start_codons_search_fwd{$chr}{$fivePrime});
                    for (@over4){
                        my ($gene,$meta_pos)=$_=~/(^[^;]+);(-?\w+$)/;
                        $start_region_signal_5{$gene}++;
                        $start_codon_meta_positions_5{$gene}{$meta_pos}+=1;
                        $start_codon_meta_positions_by_length_5{$gene}{$meta_pos}{length($seq)}+=1;
                    }
                }
            }
        }
    }                    
}
close (BAM);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#rank the genes by start and stop window expression, to exclude the bottom 10%
my %black_list_start_5;
my %black_list_start_3;

my $count=0;
my $count_exclude=0;
my $number_of_genes=keys %start_region_signal_5;
my $ten_percent=$number_of_genes*$PERCENTAGE_TO_EXCLUDE;
for my $gene ( sort { $start_region_signal_5{$a} <=> $start_region_signal_5{$b} } keys(%start_region_signal_5) ){
    if ($count <= $ten_percent){
        $black_list_start_5{$gene}=1;
        $count_exclude++;
    }
    $count++;
}


$count=0;
$count_exclude=0;
$number_of_genes=keys %start_region_signal_3;
$ten_percent=$number_of_genes*$PERCENTAGE_TO_EXCLUDE;
for my $gene ( sort { $start_region_signal_3{$a} <=> $start_region_signal_3{$b} } keys(%start_region_signal_3) ){
    if ($count <= $ten_percent){
        $black_list_start_3{$gene}=1;
        $count_exclude++;
    }
    $count++;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#scale by the counts within the start or stop window by the total counts in the window, per gene

my %barchart_upstream_5;
my %barchart_upstream_3;

for my $gene_id (keys %start_codon_meta_positions_5){
    unless (exists ($black_list_start_5{$gene_id})){  #blacklisted then remove
        for my $pos (sort { $a <=> $b } keys %{ $start_codon_meta_positions_5{$gene_id}} ){
            my $gene_sum=$start_region_signal_5{$gene_id};  
            $barchart_upstream_5{$pos}+=eval{ $start_codon_meta_positions_5{$gene_id}{$pos}/$gene_sum } || 0; 
            #barchart positions are scaled by window counts per gene
        }
    }
}

for my $gene_id (keys %start_codon_meta_positions_3){
    unless (exists ($black_list_start_3{$gene_id})){  #blacklisted then remove
        for my $pos (sort { $a <=> $b } keys %{ $start_codon_meta_positions_3{$gene_id}} ){
            my $gene_sum=$start_region_signal_3{$gene_id};
            $barchart_upstream_3{$pos}+=eval{ $start_codon_meta_positions_3{$gene_id}{$pos}/$gene_sum } || 0;
        }
    }
}

#scale lengths
my %scaled_upstream_3;
my %scaled_upstream_5;

for my $gene (sort keys %start_codon_meta_positions_by_length_5){
    unless (exists ($black_list_start_5{$gene})){  #blacklisted then remove
        for my $pos (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_5{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_5{$gene}{$pos}}){
                my $scaled_count=eval { $start_codon_meta_positions_by_length_5{$gene}{$pos}{$length}/$start_region_signal_5{$gene}} || 0 ; #divide count by the window sum for this gene
                $scaled_upstream_5{$pos}{$length}+=$scaled_count;
            }
        }
    }
}

for my $gene (sort keys %start_codon_meta_positions_by_length_3){
    unless (exists ($black_list_start_3{$gene})){  #blacklisted then remove
        for my $pos (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_3{$gene}}){
            for my $length (sort {$a <=> $b} keys %{$start_codon_meta_positions_by_length_3{$gene}{$pos}}){
                my $scaled_count=eval { $start_codon_meta_positions_by_length_3{$gene}{$pos}{$length}/$start_region_signal_3{$gene}} || 0 ; #divide count by the window sum for this gene
                $scaled_upstream_3{$pos}{$length}+=$scaled_count;
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $outStart=$outDir."/".$prefix."_start_3prime_scale.csv";
my $outStartLengthsScale=$outDir."/".$prefix."_start_lengths_scale_3prime.csv";

open (OUT1,">$outStart") || die "can't open $outStart\n";
open (OUT3,">$outStartLengthsScale") || die "can't open $outStartLengthsScale\n";

# meta plots start codons
for my $pos (sort {$a <=> $b} keys %barchart_upstream_3){
    print OUT1 "$pos,$barchart_upstream_3{$pos}\n";
}

# output scaled length values
for my $pos (sort {$a <=> $b} keys %scaled_upstream_3){ for my $length (sort {$a <=> $b} keys %{$scaled_upstream_3{$pos}}){ print OUT3 "$pos\t$length\t$scaled_upstream_3{$pos}{$length}\n"; } }

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
$outStart=$outDir."/".$prefix."_start_5prime_scale.csv";
$outStartLengthsScale=$outDir."/".$prefix."_start_lengths_scale_5prime.csv";

open (OUT5,">$outStart") || die "can't open $outStart\n";
open (OUT7,">$outStartLengthsScale") || die "can't open $outStartLengthsScale\n";

# meta plots start codons
for my $pos (sort {$a <=> $b} keys %barchart_upstream_5){
    print OUT5 "$pos,$barchart_upstream_5{$pos}\n";
}

# output scaled length values
for my $pos (sort {$a <=> $b} keys %scaled_upstream_5){ for my $length (sort {$a <=> $b} keys %{$scaled_upstream_5{$pos}}){ print OUT7 "$pos\t$length\t$scaled_upstream_5{$pos}{$length}\n"; } }

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

exit;
