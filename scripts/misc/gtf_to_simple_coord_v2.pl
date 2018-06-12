#!/usr/bin/perl -w
use strict;

#script to convert the output of Elivis repartiton perdiction to bed format.
my $inFILE=$ARGV[0];

my $out1="/export/valenfs/projects/adam/final_results/scripts/git/data_files/benchmark/tico/simple.Chromosome.coord";
my $out2="/export/valenfs/projects/adam/final_results/scripts/git/data_files/benchmark/tico/simple.pCol1B9_SL1344.coord";
my $out3="/export/valenfs/projects/adam/final_results/scripts/git/data_files/benchmark/tico/simple.pRSF1010_SL1344.coord";
my $out4="/export/valenfs/projects/adam/final_results/scripts/git/data_files/benchmark/tico/simple.pSLT_SL1344.coord";

open (OUT1,">$out1") || die;
open (OUT2,">$out2") || die;
open (OUT3,">$out3") || die;
open (OUT4,">$out4") || die;


my $count=1;

open (IN,$inFILE) || die;
while (<IN>){
    unless (/^#/){
        my @b=split("\t");
        my $chr=$b[0];
        my $class=$b[2];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];

        #my ($biotype) = $b[8] =~ /gene_biotype\s"([^\"]+)"/;

        if ($class eq "gene"){   
            my ($gene_id) = $b[8] =~ /gene_id\s"(\w+)"/;
            $gene_id=~s/_//g;

            if ($chr eq "pRSF1010_SL1344"){
                print OUT3 ">".$count."_".$start."_".$end."_".$dir."\n";
            }elsif ($chr eq "pCol1B9_SL1344"){
                print OUT2 ">".$count."_".$start."_".$end."_".$dir."\n";
            }elsif ($chr eq "pSLT_SL1344"){
                print OUT4 ">".$count."_".$start."_".$end."_".$dir."\n";
            }elsif ($chr eq "Chromosome"){
                print OUT1 ">".$count."_".$start."_".$end."_".$dir."\n";
            }
        }

        $count++;

    }
}
