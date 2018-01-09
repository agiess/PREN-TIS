#!/usr/bin/perl -w
use strict;

#script to convert the output of the tico gff of bed format.
my $refGTF=$ARGV[0];
my $predTICO=$ARGV[1];

my %gene_fwd; #chr, stop = start
my %gene_rev;

open (IN1,$refGTF) || die;
while (<IN1>){
    unless (/^#/){
        my @b=split("\t");
        my $class=$b[2];
        my $chr=$b[0];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];
        #my $gene_id="unknown";
        #($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;
        
        if ($class eq "CDS"){
        
            if ($dir eq "+"){ #fwd cases
                $gene_fwd{$chr}{$end}=$start;
            }else{ #rev cases
               $gene_rev{$chr}{$start}=$end;
           }
        }
    }
}
close(IN1);

print "track_name=tico_annotated_TIS itemRgb=On\n";

open (IN2,$predTICO) || die;
while (<IN2>){

    unless (/^#/){
        my @b=split("\t");
        my @item1=split(":",$b[0]);
        my $chr=$item1[3];
        my $start=$b[3];
        my $stop=$b[4];
        my $dir=$b[6];
        my $type="Annotated";

        if ($dir eq "+"){

            if (exists ($gene_fwd{$chr}{$stop})){
                if ($start < $gene_fwd{$chr}{$stop}){
                    $type="Extension";
                }elsif($start > $gene_fwd{$chr}{$stop}){
                    $type="Truncation";
                }
                $start-=1;
                #$stop+=3;
                print "$chr\t$start\t$stop\t$type\t1\t+\t$start\t$stop\t0,204,0\n";
            }

        }else{

            if (exists ($gene_rev{$chr}{$start})){
                if ($stop > $gene_rev{$chr}{$start}){
                    $type="Extension";
                }elsif($stop < $gene_rev{$chr}{$start}){
                    $type="Truncation";
                }
                $start-=1;
                #$stop-=1;
                print "$chr\t$start\t$stop\t$type\t1\t-\t$start\t$stop\t0,204,0\n";
            }

        }
    }      
}
close(IN2);

exit;
