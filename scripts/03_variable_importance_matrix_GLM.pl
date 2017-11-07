#!/usr/bin/perl -w
use strict;

#30/11/17
my $inFile=$ARGV[0];
my $outFile=$ARGV[1];

my $upstream_TIS=0;   #canonical_candidate_sites_ORF_region
my $reads_up=0;       #proprtion_of_reads_upstream
my $reads_down=0;     #proportion_of_reads_downstrea
my $reads_position=0; #proportion_of_reads_at_positionm
my $ORF_FPKM=0;       #ORF_fpkm

my %length_matrix;
for my $pos (-20 .. 20){
   for my $length (20 .. 49){
       $length_matrix{$pos}{$length}=0;
   }
}

my $best_codon=0;
my $best_codon_positive=0;
my %best_seq;
my %best_seq_positive;
for (-20 .. 20){
   $best_seq{$_}=0;
   $best_seq_positive{$_}=0;
}

open (IN,$inFile) || die "can't open $inFile\n";
while (<IN>){

    unless (/^"",/){ #skip header

        chomp();
        my @line=split(",");

        my $name=$line[1];
        my $coefficient=$line[2];
        my $scaled_coefficient=$line[3];

        #strip quotes
        $name =~ tr/"/'/d;   
        $coefficient =~ tr/"/'/d;
        $scaled_coefficient =~ tr/"/'/d;
 
        #names
        if ($name =~ /canonical_candidate_sites_ORF_region/ ){ $upstream_TIS=$scaled_coefficient; }
        if ($name =~ /proportion_of_reads_upstream/ ){ $reads_up=$scaled_coefficient; }                  
        if ($name =~ /proportion_of_reads_downstream/ ){ $reads_down=$scaled_coefficient; }
        if ($name =~ /proportion_of_reads_at_position/ ){ $reads_position=$scaled_coefficient; }
        if ($name =~ /window_up_down_ratio/ ){ $reads_position=$scaled_coefficient; }                           #to remove
        if ($name =~ /ORF_fpkm/ ){ $ORF_FPKM=$scaled_coefficient; }

        elsif ($name =~ /codon\.\w{3}/ ){ 
            my $tmp_coef=$scaled_coefficient; #check for -ve/+ve
            if ($scaled_coefficient < 0){ $tmp_coef=$scaled_coefficient*-1; }    
            if ($tmp_coef > $best_codon_positive ){
                $best_codon=$scaled_coefficient;
                $best_codon_positive=$tmp_coef;
            }
        }

        elsif ($name =~ /seq_(\w+)\.\w/ ){
            my $tmp_coef=$scaled_coefficient; #check for -ve/+ve
            if ($scaled_coefficient < 0){ $tmp_coef=$scaled_coefficient*-1; }              
            if ($tmp_coef > $best_seq{$1} ){
                $best_seq{$1}=$scaled_coefficient;
                $best_seq_positive{$1}=$tmp_coef;
            }
        }

        elsif ($name =~ /seq_\.(\w+)\.\w/ ){
            my $pos=-$1;
            my $tmp_coef=$scaled_coefficient; #check for -ve/+ve
            if ($scaled_coefficient < 0){ $tmp_coef=$scaled_coefficient*-1; }
            if ($tmp_coef > $best_seq{$pos} ){
                $best_seq{$pos}=$scaled_coefficient;
                $best_seq_positive{$pos}=$tmp_coef;
            }
        }

        elsif ( $name =~ /X(\d+)\_(\d+)/ ){
            $length_matrix{$1}{$2}=$scaled_coefficient;
        }

        elsif ( $name =~ /X\.(\d+)\_(\d+)/ ){
            my $pos=-$1;
            $length_matrix{$pos}{$2}=$scaled_coefficient;
        } 
    }	
}
close(IN);

#add the sequence to the matrix as length 19nt
for (sort {$a <=> $b} keys %best_seq){
    $length_matrix{$_}{19}=$best_seq{$_};
}

#add the codon to the matrix as length 19nt, positions0-2nt
$length_matrix{0}{19}=$best_codon;
$length_matrix{1}{19}=$best_codon;
$length_matrix{2}{19}=$best_codon;

#add the proportion of reads in the window to the matrix as length 50nt.
for (-20 .. -1){
   $length_matrix{$_}{50}=$reads_up;
}
for (1 .. 20){
    $length_matrix{$_}{50}=$reads_down;
}
$length_matrix{0}{50}=$reads_position;

#add ORF_FMPK to matrix as length 51
for (-20 .. 20){
   $length_matrix{$_}{51}=$ORF_FPKM;
}

#add count of upstream TIS to matrix as length 52
for (-20 .. 20){
   $length_matrix{$_}{52}=$upstream_TIS;
}

open (OUT, ">$outFile") || die "can't open $outFile\n";

#print in long format
print OUT  "position,length,standardized_coefficients\n";
for my $pos (sort {$a <=> $b} keys %length_matrix){
    for my $length (sort {$a <=> $b} keys %{ $length_matrix{$pos} }){
        print OUT "$pos,$length,$length_matrix{$pos}{$length}\n";
    }
}

exit;
