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
my %best_seq;
for (-20 .. 20){
   $best_seq{$_}=0;
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
        $name =~ s/\"//g;
        $coefficient =~ s/\"//g;
        $scaled_coefficient =~ s/\"//g;
       
        #names
        if ($name =~ /distance_to_region_start/ ){ $upstream_TIS=$scaled_coefficient; }
        if ($name =~ /proportion_of_reads_upstream/ ){ $reads_up=$scaled_coefficient; }                  
        if ($name =~ /proportion_of_reads_downstream/ ){ $reads_down=$scaled_coefficient; }
        if ($name =~ /proportion_of_reads_at_position/ ){ $reads_position=$scaled_coefficient; }
        if ($name =~ /ORF_fpkm/ ){ $ORF_FPKM=$scaled_coefficient; }

        elsif ($name =~ /codon/ ){  
            if ( $scaled_coefficient > $best_codon ){
                $best_codon=$scaled_coefficient;
            }
        }

        elsif ($name =~ /seq_(\w+)/ ){
            if ( $scaled_coefficient > $best_seq{$1} ){
                 $best_seq{$1}=$scaled_coefficient;
            }
        }

        elsif ($name =~ /seq_\.(\w+)/ ){
            my $pos=-$1;
            if ( $scaled_coefficient > $best_seq{$pos} ){
                $best_seq{$pos}=$scaled_coefficient;
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
print OUT "position,length,value\n"; 
for my $pos (sort {$a <=> $b} keys %length_matrix){
    for my $length (sort {$a <=> $b} keys %{ $length_matrix{$pos} }){
        print OUT "$pos,$length,$length_matrix{$pos}{$length}\n";
    }
}

exit;
