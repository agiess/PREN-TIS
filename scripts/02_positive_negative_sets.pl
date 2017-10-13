#!/usr/bin/perl -w
use strict;

#to do 29/09/2017
#script to take a matrix of genomic positions and an optional bed file of supported start sites (for example N-terminal proteomics)
#assign codons to neative and positive sets
#take only the 50% highest genes by cds expression in thepositive set

my $matrix=$ARGV[0];
my $out_positive=$ARGV[1]; #1000 start codons from upper 50% of genes by experssions, without n-termial support
my $out_negative=$ARGV[2]; #8000 non start codons, without n-termini support (inframe)
my $bed_file=$ARGV[3];     #hi confidence start sites (optional)

###################################################################################################
my %n_term_fwd; #key=chr, key=start_position, value=count
my %n_term_rev;

if ($bed_file){
    #open bed file

    open (BED, $bed_file) || die "can't open $bed_file";      #bed is O based at start, 1 based at end
    while (<BED>){

        unless (/^track/){         #skip header
            my @b=split("\t");
	     	my $chr=$b[0];
            my $start=$b[1]+1;     #start is zero bases
            my $stop=$b[2];
            my $count=$b[4];
            my $dir=$b[5];
                    
            #assign to metaplots 
            if ($dir eq "+"){                         #fwd cases
                #five_prime=$start;
                $n_term_fwd{$chr}{$start}+=$count;            
            }else{                                    #reverse cases
                #five_prime=$stop;
                 $n_term_rev{$chr}{$stop}+=$count;
            }
		}
    }                    
}
close (BED);

###################################################################################################
#open matrix

my @start_supported;
my @start_random;
my @start_FPKM; #for median
my @not_start;
my $header="#empty_header\n";

open (MAT, $matrix) || die "can't open $matrix";     
while (<MAT>){

    if (/^#/){
		$header=$_;
	}else{

        #1  id
        #2  codon
        #3  dir
        #4  codon_rank
        #5  annotated_start_site
        #6  near_cognate
        #7  reads_at_pos
        #8  window_reads_downstream
        #9  window_reads_upstream
        #10 window_up_down_ratio
        #11 proportion_downstream
        #12 proportion_upstream
        #13 ORF_FPKM

        my @l=split(",");
        my $dir=$l[2];
        my $start_codon=$l[4];
        my $fpkm=$l[12]; 

		my ($chr, $pos) = $l[0] =~ /^(.*)_(\d+)_(fwd|rev)$/;     

        #start codons
        if ($start_codon eq "TRUE"){ 
    
            if ($dir eq "fwd"){

                if (exists ($n_term_fwd{$chr}{$pos})){
                    push (@start_supported,$_);  #for counting
                    push (@start_FPKM,$fpkm); 
                }else{
                    push (@start_random,$_);
                    push (@start_FPKM,$fpkm);
                }

            }else{ #rev
                if (exists ($n_term_rev{$chr}{$pos})){
                    push (@start_supported,$_);   #for counting
                    push (@start_FPKM,$fpkm);
                }else{
                    push (@start_random,$_);
                    push (@start_FPKM,$fpkm);
                }
            }

        }else{ #non start codons
            if ($dir eq "fwd"){
                unless (exists ($n_term_fwd{$chr}{$pos})){
                    push (@not_start,$_);
                }
            }else{ #rev
                unless (exists ($n_term_rev{$chr}{$pos})){
                    push (@not_start,$_);
                }
            }
    	}
    }
}
close(MAT);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#find median expression from start codons
my $median_exp=0;
my $size=$#start_FPKM;
my $count=0;
for my $ORF (sort {$a <=> $b} @start_FPKM){
	$count++;
	if ($count > $size/2){
        $median_exp=$ORF;
		last;
	}
}

print "size = $size, median = $median_exp\n";
my $sizeUn=$#start_random; 
my $sizeSu=$#start_supported; 

if ($bed_file){
    if ($sizeUn < 1800){
        print "ERROR: There are not enough unsupported start codons to train the model, the predictions script requires at least 1700 positive examples\n";
        exit;
    }else{
        print "there are $sizeUn unsupported start codons\n";
        print "there are $sizeSu supporterd start codons\n";
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#select all start codons from the upper 50% of genes by expression, without N-termini support
my @positive;
for my $line (@start_random){
    my @position=split(",",$line);
    if ($position[12] >= $median_exp){ #filter on median exp!!!!!
    	push @positive, $line;
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output negative set3
my @negative;
my $size2=0;
while ($size2 < 8000){
 
    #get element
    my $index = rand @not_start;
    push (@negative, splice @not_start, $index, 1);
    $size2++;
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output
open (OUT1,">$out_positive")  || die "can't open $out_positive\n";
open (OUT2,">$out_negative")  || die "can't open $out_negative\n";

#headers
print OUT1 $header;
print OUT2 $header;

for (@positive){
    print OUT1 $_;
}

for (@negative){
    print OUT2 $_;
}

exit;