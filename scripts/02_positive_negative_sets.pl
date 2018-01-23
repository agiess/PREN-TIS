#!/usr/bin/perl -w
use strict;

#07/12/2017
#script to take a matrix of genomic positions and an optional bed file of supported start sites (for example N-terminal proteomics)
#assign codons to neative and positive sets
#take the 50% highest genes by cds expression in the positive set (exluding positions that will later be used for validation)
#take a randonly selected negative set of 4 times the number of positive examples (exluding poitions that will later be used for validation)

my $matrix=$ARGV[0];
my $out_train_positive=$ARGV[1]; #Start codons from upper 50% of genes by expression,
my $out_train_negative=$ARGV[2]; #Inframe non-start codons
my $out_test_positive=$ARGV[3];  #Start codons from upper 50% of genes by expression,
my $out_test_negative=$ARGV[4];  #Inframe non-start codons
my $proportion=$ARGV[5];         #the proportion of the 50% highest genes to select for positives
my $bed_file=$ARGV[6];           #validated start sites (optional)

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
                $n_term_fwd{$chr}{$start}=1;            
            }else{                                    #reverse cases
                #five_prime=$stop;
                 $n_term_rev{$chr}{$stop}=1;
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
        #my $start_codon=$l[4];
        #my $fpkm=$l[12]; 
        my $start_codon=$l[13];
        my $fpkm=$l[21];

        my ($chr, $pos) = $l[0] =~ /^(.*)_(\d+)_(fwd|rev)/;
 
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
my $size=@start_FPKM;
my $count=0;
for my $ORF (sort {$a <=> $b} @start_FPKM){
	$count++;
	if ($count > $size/2){
        $median_exp=$ORF;
		last;
	}
}

print "There are $size genes\n";
print "The median expression value is $median_exp FPKM\n";
print "The proportion of highly expressed genes to select is $proportion\n";

my $positive_train=int((($size/2)*$proportion)*0.8);
my $positive_test=int((($size/2)*$proportion)*0.2);

if ($bed_file){
    my $sizeUn=@start_random;
    my $sizeSu=@start_supported;

    print "There are $sizeUn start codons that are not supported by the validation set\n";
    print "There are $sizeSu start codons supported by the validation set\n";


    #check if there are enough highly expressed start codons after exlcuding validated positions
    if ((($size/2)-$sizeSu) < ($positive_train+$positive_test) ){
        print "there are not enough positive examples after exluding validated positions, either use a smaller validation set or lower proportion of highly expressed genes\n";
        exit 1;
    }
}

#set seed
srand(7777);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#select all start codons from the upper 50% of genes by expression, without N-termini support
my @positive_train;
my @positive_test;
my $size1=0;

while ($size1 < ($positive_train+$positive_test)){

    my $index = rand @start_random;
    my @position=split(",",$start_random[$index]);
#    if ($position[12] >= $median_exp){ #filter on median expression
    if ($position[21] >= $median_exp){ #filter on median expression

        if ($size1 < $positive_train){
            push (@positive_train, splice @start_random, $index, 1);
        }else{
            push (@positive_test, splice @start_random, $index, 1);
        }
        $size1++;
    }
}

my $pos_train_size=@positive_train;
my $pos_test_size=@positive_test;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#The negative sets should have 4 times the number of positive examples
my @negative_train;
my @negative_test;
my $size2=0;
while ($size2 < (($positive_train+$positive_test)*4)){
 
    my $index = rand @not_start;
    if ($size2 < ($positive_train*4)){
       push (@negative_train, splice @not_start, $index, 1);
    }else{
       push (@negative_test, splice @not_start, $index, 1);
    }
    $size2++;
}

my $neg_train_size=@negative_train;
my $neg_test_size=@negative_test;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#output
open (OUT1,">$out_train_positive")  || die "can't open $out_train_positive\n";
open (OUT2,">$out_train_negative")  || die "can't open $out_train_negative\n";
open (OUT3,">$out_test_positive")  || die "can't open $out_test_positive\n";
open (OUT4,">$out_test_negative")  || die "can't open $out_test_negative\n";

#headers
print OUT1 $header;
print OUT2 $header;
print OUT3 $header;
print OUT4 $header;

for (@positive_train){
    print OUT1 $_;
}

for (@negative_train){
    print OUT2 $_;
}

for (@positive_test){
    print OUT3 $_;
}

for (@negative_test){
    print OUT4 $_;
}

print "There are $pos_train_size positions in the positive training set\n";
print "There are $neg_train_size positions in the negative training set\n";
print "There are $pos_test_size positions in the positive testing set\n";
print "There are $neg_test_size positions in the negative testing set\n";

exit;
