#!/usr/bin/perl -w
#use strict;

#script to compare the output of the triTISA med2.0 file to reference and convert to bed format.
my $refGTF=$ARGV[0];
my $predMED=$ARGV[1];

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

#chr from name
my $chr="Chromosome";
if ($predMED =~/pRSF1010/){ $chr="pRSF1010_SL1344"; }
if ($predMED =~/pCol1B9/){ $chr="pCol1B9_SL1344"; }
if ($predMED =~/pSLT/){ $chr="pSLT_SL1344"; }
if ($predMED =~/Chromosome/){ $chr="Chromosome"; }

print "track_name=elvis_annotated_TIS itemRgb=On\n";

open (IN2,$predMED) || die;
while (<IN2>){

    unless (/^#/){
#       my @b=split("\t");
        #my $chr="Chromosome"; #can I parse this fron the 
        my ($start,$stop,$dir)=$_=~/^\s+([^\s]+)\s+([^\s]+)\s([^\s]+)/;

#       my $start=$b[0];
#       my $stop=$b[1];
#       my $dir=$b[2];
        my $type="Annotated";

#        print "$start,$stop,$dir\n";

        if ($dir eq "+"){
           
            if (exists ($gene_fwd{$chr}{$stop})){
                if ($start < $gene_fwd{$chr}{$stop}){
                    $type="Extension";
                }elsif($start > $gene_fwd{$chr}{$stop}){
                    $type="Truncation";
                }
                $start-=1;
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
                print "$chr\t$start\t$stop\t$type\t1\t-\t$start\t$stop\t0,204,0\n";
            }
        }
    }      
}
close(IN2);

exit;
