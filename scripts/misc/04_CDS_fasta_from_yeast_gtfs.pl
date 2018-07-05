#!/usr/bin/perl -w
use strict;

#to do 03/07/2018
#script to produce CDS fasta regions from GTF + updated predictions (updated to include elongated TIS)

#input files:
my $gtf=$ARGV[0];
my $fasta=$ARGV[1];
my $pred_bed=$ARGV[2];

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open gtf and get transcript lengths
my %transcripts; #key = gene_id, transcript_id, #value = sum_exon_lengths;

my %all_annotated_start_sites; #for checking if predictions are novel TIS or isoform varients

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
            if ($class eq "start_codon"){
                if ($dir eq "+"){ #fwd cases. Use start positions as 5'
                    $all_annotated_start_sites{$start}=1;
                }else{
                    $all_annotated_start_sites{$end}=1;
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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#parse bed and save elongations here
my %extended_genes_fwd;
my %extended_genes_rev;
my %extended_gene_start_pos_fwd;
my %extended_gene_start_pos_rev;
my %extended_gene_start_coord_fwd;
my %extended_gene_start_coord_rev;

open (BED, $pred_bed) || die;
while (<BED>){
    unless (/^track/){
        my @line=split("\t");
        my $chr=$line[0];
        my $start=$line[1]+1;
        my $stop=$line[2];
        my $class=$line[3];
        my $dir=$line[5];

        if ($class eq "Extension"){
            if ($dir eq "+"){
                $extended_genes_fwd{$chr}{$stop}=$start;
            }else{
                $extended_genes_rev{$chr}{$start}=$stop;
            }
        }
    }
}
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#setup transcript models

my %gene_model_fwd;
my %start_coord_fwd;
my %stop_coord_fwd;
my %three_prime_most_coord_fwd;
my %five_prime_most_coord_fwd;

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon
        my $chr= $gene_2_chr{$gene};
        my $model_pos=0;
        $five_prime_most_coord_fwd{$gene}=$model_pos;  #initalise the 5' to the first coord coord
        for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
            my $exon_end=$gene_exons_fwd{$gene}{$exon_start};
            for ($exon_start .. $exon_end){
                $gene_model_fwd{$gene}{$model_pos}=$_;
                if ($_ == $gene_stop_codon_fwd{$gene}){
                    $stop_coord_fwd{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($_ == $gene_start_codon_fwd{$gene}){
                    $start_coord_fwd{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                if (exists ($extended_genes_fwd{$chr}{$_})){
                    $extended_gene_start_pos_fwd{$gene}=$extended_genes_fwd{$chr}{$_};
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
my %three_prime_most_coord_rev;
my %five_prime_most_coord_rev;

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon
        my $chr= $gene_2_chr{$gene};
        my $model_pos=0;
        $five_prime_most_coord_rev{$gene}=$model_pos;  #initalise the 5' to the first coord coord
        for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
            my $exon_start=$gene_exons_rev{$gene}{$exon_end};
            while ($exon_start >= $exon_end){
                $gene_model_rev{$gene}{$model_pos}=$exon_start;
                if ($exon_start == $gene_stop_codon_rev{$gene}){
                    $stop_coord_rev{$gene}=$model_pos;    #find the index of the stop codon per gene
                }
                if ($exon_start == $gene_start_codon_rev{$gene}){
                    $start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
                }
                if (exists ($extended_genes_rev{$chr}{$exon_start})){
                    $extended_gene_start_pos_rev{$gene}=$extended_genes_rev{$chr}{$exon_start};
                }
                $model_pos++;
                $exon_start--;
            }
        }
        $three_prime_most_coord_rev{$gene}=$model_pos-1; #store the 3 prime most position of each gene
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

for my $gene (keys %gene_exons_fwd){
    if ( (exists ($gene_start_codon_fwd{$gene})) && (exists ($gene_stop_codon_fwd{$gene})) ) { #restrict to genes with annotated start + stop codon

       if (exists ($extended_gene_start_pos_fwd{$gene})){

            my $model_pos=0;
            for my $exon_start (sort {$a <=> $b} keys %{ $gene_exons_fwd{$gene} } ){
                my $exon_end=$gene_exons_fwd{$gene}{$exon_start};
                for ($exon_start .. $exon_end){

                    if ($_ == $extended_gene_start_pos_fwd{$gene}){
                        $extended_gene_start_coord_fwd{$gene}=$model_pos;
#                       print "$gene,fwd,extended,$_,$model_pos\n";
                    }
                    $model_pos++;
                }
            }
        }
    }
}

for my $gene (keys %gene_exons_rev){
    if ( (exists ($gene_start_codon_rev{$gene})) && (exists ($gene_stop_codon_rev{$gene})) ) { #restrict to genes with annotated start + stop codon
        if (exists ($extended_gene_start_pos_rev{$gene})){
            my $model_pos=0;
            for my $exon_end (reverse (sort {$a <=> $b} keys %{ $gene_exons_rev{$gene} } )){
                my $exon_start=$gene_exons_rev{$gene}{$exon_end};
                while ($exon_start >= $exon_end){

                    if ($exon_start == $extended_gene_start_pos_rev{$gene}){
                        $extended_gene_start_coord_rev{$gene}=$model_pos;    #find the index of the start codon per gene
#                       print "$gene,rev,extended,$exon_start,$model_pos\n";
                    }
                    $model_pos++;
                    $exon_start--;
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

for my $gene (keys %gene_model_fwd){
    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_fwd{$gene};
    my $stop_coord=$stop_coord_fwd{$gene};

    if (exists ($extended_gene_start_coord_fwd{$gene})){
        $start_coord=$extended_gene_start_coord_fwd{$gene}; #update extended start codon coord
    }

    my @CDS;
    for my $coord (sort {$a <=> $b} keys %{ $gene_model_fwd{$gene} } ){
        if (($coord >= ($start_coord)) && ($coord <= ($stop_coord-1))){ #get the CDS (exluding stop sodon)
            if (exists ($gene_model_fwd{$gene}{$coord})){
                push (@CDS, substr($fasta_sequences{$chr}, ($gene_model_fwd{$gene}{$coord} -1), 1) );
            }else{
                push (@CDS, "N");
            }
        }
    }
    my $nuc_CDS=join("", @CDS);
    my $pro_CDS=&translate($nuc_CDS);
    print ">$gene\n$pro_CDS\n";
}

for my $gene (keys %gene_model_rev){

    my $chr=$gene_2_chr{$gene};
    my $start_coord=$start_coord_rev{$gene};
    my $stop_coord=$stop_coord_rev{$gene};

    if (exists ($extended_gene_start_coord_rev{$gene})){
        $start_coord=$extended_gene_start_coord_rev{$gene}; #update extended start codon coord
    }

    my @CDS;
    for my $coord (sort {$a <=> $b} keys %{ $gene_model_rev{$gene} } ){
        if (($coord >= ($start_coord)) && ($coord <= ($stop_coord-1))){ #get the CDS (excluding stop codon)
            if (exists ($gene_model_rev{$gene}{$coord})){
                push (@CDS, substr($fasta_sequences{$chr}, ($gene_model_rev{$gene}{$coord} -1), 1) );
            }else{
                push (@CDS, "N");
            }
        }
    }
    my $nuc_CDS=join("", @CDS);
    $nuc_CDS=~tr/ACGTacgt/TGCAtgca/;
    my $pro_CDS=&translate($nuc_CDS);
    print ">$gene\n$pro_CDS\n";
}

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub translate{
    #translate a the nuncleotides of a CDS into its corrosponding amino acid sequence
 
    my $nucleotice_sequence=$_[0];

    my(%genetic_code) = (
      'TCA' => 'S',    # Serine
      'TCC' => 'S',    # Serine
      'TCG' => 'S',    # Serine
      'TCT' => 'S',    # Serine
      'TTC' => 'F',    # Phenylalanine
      'TTT' => 'F',    # Phenylalanine
      'TTA' => 'L',    # Leucine
      'TTG' => 'L',    # Leucine
      'TAC' => 'Y',    # Tyrosine
      'TAT' => 'Y',    # Tyrosine
      'TAA' => '_',    # Stop
      'TAG' => '_',    # Stop
      'TGC' => 'C',    # Cysteine
      'TGT' => 'C',    # Cysteine
      'TGA' => '_',    # Stop
      'TGG' => 'W',    # Tryptophan
      'CTA' => 'L',    # Leucine
      'CTC' => 'L',    # Leucine
      'CTG' => 'L',    # Leucine
      'CTT' => 'L',    # Leucine
      'CCA' => 'P',    # Proline
      'CCC' => 'P',    # Proline
      'CCG' => 'P',    # Proline
      'CCT' => 'P',    # Proline
      'CAC' => 'H',    # Histidine
      'CAT' => 'H',    # Histidine
      'CAA' => 'Q',    # Glutamine
      'CAG' => 'Q',    # Glutamine
      'CGA' => 'R',    # Arginine
      'CGC' => 'R',    # Arginine
      'CGG' => 'R',    # Arginine
      'CGT' => 'R',    # Arginine
      'ATA' => 'I',    # Isoleucine
      'ATC' => 'I',    # Isoleucine
      'ATT' => 'I',    # Isoleucine
      'ATG' => 'M',    # Methionine
      'ACA' => 'T',    # Threonine
      'ACC' => 'T',    # Threonine
      'ACG' => 'T',    # Threonine
      'ACT' => 'T',    # Threonine
      'AAC' => 'N',    # Asparagine
      'AAT' => 'N',    # Asparagine
      'AAA' => 'K',    # Lysine
      'AAG' => 'K',    # Lysine
      'AGC' => 'S',    # Serine
      'AGT' => 'S',    # Serine
      'AGA' => 'R',    # Arginine
      'AGG' => 'R',    # Arginine
      'GTA' => 'V',    # Valine
      'GTC' => 'V',    # Valine
      'GTG' => 'V',    # Valine
      'GTT' => 'V',    # Valine
      'GCA' => 'A',    # Alanine
      'GCC' => 'A',    # Alanine
      'GCG' => 'A',    # Alanine
      'GCT' => 'A',    # Alanine
      'GAC' => 'D',    # Aspartic Acid
      'GAT' => 'D',    # Aspartic Acid
      'GAA' => 'E',    # Glutamic Acid
      'GAG' => 'E',    # Glutamic Acid
      'GGA' => 'G',    # Glycine
      'GGC' => 'G',    # Glycine
      'GGG' => 'G',    # Glycine
      'GGT' => 'G',    # Glycine
    );

    my $amino_acids;
    my @codons = $nucleotice_sequence =~ /(.{1,3})/g;
 
    #my @codons;
    #push @codons, substr $nucleotice_sequence, 0, 3, '' while $nucleotice_sequence;

    for (@codons){
        if(exists $genetic_code{$_}){
            $amino_acids.=$genetic_code{$_};
        }else{
            print STDERR "Bad codon \"$_\"!!\n";
        }
    }
    return $amino_acids;
}
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

