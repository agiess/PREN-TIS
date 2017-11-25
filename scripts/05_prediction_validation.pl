#!/usr/bin/perl -w
use strict;

#to do 20/10/2017
#script to take a bed file of orf predictions and to report n-termini suppoprt per prediction catagory
#looking for exact matches only

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $nterm=$ARGV[2]; #bed of N-termini predictions
my $bed=$ARGV[3];   #of predictions

#for all matching, elongated or predicted
#count the number that have n-terminal support as matching, elongated or predicted

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open gtf file setup genomic hashes, to store start and stop coords per gene, per direction
my %gene_regions_fwd_start; #old way of finding starts and ends
my %gene_regions_fwd_stop;
my %gene_regions_rev_start;
my %gene_regions_rev_stop;

#get gene coords per direction
my %start_codon; #key=gene_id, value=pos
my %gene_regions_fwd; # chr, start, stop = gene 
my %gene_regions_rev; # chr, start, stop = gene 

#store gene ids
my %gene_chr;  #gene = chr
my %gene_fwd_start; #gene = start
my %gene_fwd_stop ; #gene = end 
my %gene_rev_start; #gene = start 
my %gene_rev_stop ; #gene = end

#for subroutine
my %gene_info_fwd; #key=gene_id, key2=start, value=stop
my %gene_info_rev; #key=gene_id, key2=start, value=stop

open(GENES,$gtf) || die "can't open $gtf";        #gtf is 1 based
while (<GENES>){
    unless(/^#/){
        my @b=split("\t");
        my $class=$b[2];
        my $chr=$b[0];
        my $prim5=$b[3];
        my $prim3=$b[4];
        my $dir=$b[6];
        my $gene_id="unknown";
        ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;

        #update to use cds rather than start / stop
        if ($class eq "CDS"){

            if ($dir eq "+"){ #use start positions as begining of feature
                $gene_regions_fwd_start{$gene_id}=$prim5;
                $gene_regions_fwd_stop{$chr}{$prim3}=$gene_id;
                $gene_regions_fwd{$chr}{$prim5}{$prim3}=$gene_id;
                $gene_fwd_start{$gene_id}=$prim5;
                $gene_fwd_stop{$gene_id}=$prim3;
                $gene_chr{$gene_id}=$chr;
                $gene_info_fwd{$gene_id}{$prim5}=$prim3;
            }else{
                $gene_regions_rev_start{$gene_id}=$prim3;
                $gene_regions_rev_stop{$chr}{$prim5}=$gene_id;
                $gene_regions_rev{$chr}{$prim5}{$prim3}=$gene_id;
                $gene_rev_start{$gene_id}=$prim5;
                $gene_rev_stop{$gene_id}=$prim3;
                $gene_chr{$gene_id}=$chr;
                $gene_info_rev{$gene_id}{$prim3}=$prim5-1;
            }
        }
    }
}
close(GENES);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open n-termini bed
my %n_term_fwd; #key=chr, key=start_position, value=count
my %n_term_rev;

open (PEP, $nterm) || die "can't open $nterm";      #bed is O based at start, 1 based at end
while (<PEP>){
    unless (/^track/){  #skip header
        my @b=split("\t");
        my $chr=$b[0];
        my $start=$b[1]+1;     #start is zero bases
        my $stop=$b[2];
        my $count=$b[4];
        my $dir=$b[5];

        #assign to metaplots 
        if ($dir eq "+"){                                      #fwd cases
            #five_prime=$start;
            $n_term_fwd{$chr}{$start}=$stop;
        }else{                                                 #reverse cases
            #five_prime=$stop;
            $n_term_rev{$chr}{$stop}=$start;
        }
    }
}
close (PEP);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open fasta 
my %fasta_sequences; #key = sequence_name, value=sequence
my $name;
open (FA, $fasta) || die "can't open $fasta";
while (<FA>){
    chomp;
    if (/^>([^\s]+)/){ #take header up to the first space
        $name=$1;
    }else{
        $fasta_sequences{$name}.=$_;
    }
}
close(FA);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#pass hash reference to subroutine
my ($inframe_fwd_ref, $inframe_rev_ref)=&stopToStopFromGTF( \%fasta_sequences, \%gene_info_fwd, \%gene_info_rev, \%gene_chr);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign N-termini here
#get a hash of genes and how they are supported
#to dereferece (makes a new copy og the hash)
my %inframe_fwd=%{$inframe_fwd_ref};  #gene,chr,pos=1
my %inframe_rev=%{$inframe_rev_ref};  #gene,chr,pos=1
my $n_sum_ano=0;
my $n_sum_tru=0;
my $n_sum_elo=0;
my $n_sum_new=0;
my %n_term_assignment_fwd; #chr, start, stop = ann/tru/elo
my %n_term_assignment_rev; #chr, stop, start = ann/tru/elo

my %supported_genes;

for my $c (sort keys %n_term_fwd ){
    for my $five (sort keys %{$n_term_fwd{$c}}){
        my $three=$n_term_fwd{$c}{$five};
        #check if they overlap a gene region and are in frame
        if (exists ($inframe_fwd{$c}{$five})){
            for my $gene (keys %{$inframe_fwd{$c}{$five}}){
                $supported_genes{$gene}=1;
                if ($five == $gene_fwd_start{$gene}){              #annotated
                    $n_term_assignment_fwd{$c}{$three}{$five}="ann";
                    $n_sum_ano++;
                }elsif($five > $gene_fwd_start{$gene} ){           #trunaction
                    $n_term_assignment_fwd{$c}{$three}{$five}="tru";
                    $n_sum_tru++;
                }elsif($five < $gene_fwd_start{$gene}){            #elongation
                    $n_term_assignment_fwd{$c}{$three}{$five}="elo";
                    $n_sum_elo++;
                }
            }
        }
    }
}

for my $c (sort keys %n_term_rev ){
    for my $five (sort keys %{$n_term_rev{$c}}){
        my $three=$n_term_rev{$c}{$five};
        #check if they overlap a gene region and are in frame
        if (exists ($inframe_rev{$c}{$five})){
            for my $gene (keys %{$inframe_rev{$c}{$five}}){
                $supported_genes{$gene}=1;
                if ($five == $gene_rev_stop{$gene}){               #annotated
                    $n_term_assignment_rev{$c}{$three}{$five}="ann";
                    $n_sum_ano++;
                }elsif($five < $gene_rev_stop{$gene} ){            #truncation
                    $n_term_assignment_rev{$c}{$three}{$five}="tru";
                    $n_sum_tru++;
                }elsif($five > $gene_rev_stop{$gene}){             #elongation                      
                    $n_term_assignment_rev{$c}{$three}{$five}="elo";
                    $n_sum_elo++;
                }
            }
        }
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#count potential TIS here.
#need a list of genes with support first
my ($total_tis)=&countTIS( \%fasta_sequences, \%gene_info_fwd, \%gene_info_rev, \%gene_chr, \%supported_genes);
my $count_supported=keys %supported_genes;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open bed of predictions
my %call; 
my $pred_sum_ano=0;
my $pred_sum_tru=0;
my $pred_sum_elo=0;

open (BED, $bed) || die;
while (<BED>){

    unless (/^track/){  #skip header
        my @b=split("\t");
        my $chr=$b[0];
        my $start=$b[1]+1;     #start is zero based
        my $stop=$b[2];
        my $count=$b[4];
        my $dir=$b[5];
        my $type=$b[3]; #Annotated Truncation Extension

        if ($type eq "Annotated"){  $pred_sum_ano++; }
        if ($type eq "Extension"){  $pred_sum_elo++; }
        if ($type eq "Truncation"){ $pred_sum_tru++; }

        if ($dir eq "+"){

            #check if the region shares a stop codon with a peptide
            if (exists ($n_term_assignment_fwd{$chr}{$stop}) ){  

                for my $nstart (keys %{ $n_term_assignment_fwd{$chr}{$stop} } ){

                    if ($type eq "Annotated"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){
                            if ($start == $nstart){ $call{$chr}{$stop}="ano_ano"; }
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){ 
                            $call{$chr}{$stop}="ano_tru"; 
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            $call{$chr}{$stop}="ano_elo";
                        }

                    #check for extensions
                    }elsif ( $type eq "Extension"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){ 
                            $call{$chr}{$stop}="elo_ano";
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){ 
                            $call{$chr}{$stop}="elo_tru"; 
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            if ($start == $nstart){ $call{$chr}{$stop}="elo_elo"; }
                        } 

                    #otherwise truncation
                    }elsif ( $type eq "Truncation"){
                        if ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "ann"){ 
                            $call{$chr}{$stop}="tru_ano";
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "tru"){
                            if ($start == $nstart){ $call{$chr}{$stop}="tru_tru"; }
                        }elsif ($n_term_assignment_fwd{$chr}{$stop}{$nstart} eq "elo"){ 
                            $call{$chr}{$stop}="tru_elo";
                        }
                    }
                }
            }
        #rev
        }else{

            #for rev, start is 3' (stop)
            #stop is 5' (start)

            #check if the region shares a stop codon with a peptide
            if (exists ($n_term_assignment_rev{$chr}{$start} ) ){

                for my $nstop (keys %{ $n_term_assignment_rev{$chr}{$start}} ){

                    if ($type eq "Annotated"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ann"){ 
                            if ($stop == $nstop){ $call{$chr}{$start}="ano_ano"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){ 
                            $call{$chr}{$start}="ano_tru";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            $call{$chr}{$start}="ano_elo";
                        }

                    }elsif ($type eq "Extension"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            if ($stop == $nstop){ $call{$chr}{$start}="elo_elo"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){ 
                            $call{$chr}{$start}="elo_tru";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ano"){ 
                            $call{$chr}{$start}="elo_ano"; 
                        }
                 
                    }elsif ($type eq "Truncation"){
                        if ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "tru"){
                            if ($stop == $nstop){ $call{$chr}{$start}="tru_tru"; }
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "ann"){ 
                            $call{$chr}{$start}="tru_ano";
                        }elsif ($n_term_assignment_rev{$chr}{$start}{$nstop} eq "elo"){ 
                            $call{$chr}{$start}="tru_elo"; 
                        }
                    }
                }
            }   
        }
    }
}
close(BED);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my $ano_ano=0;
my $ano_tru=0;
my $ano_elo=0;
my $tru_ano=0;
my $tru_tru=0;
my $tru_elo=0;
my $elo_ano=0;
my $elo_tru=0;
my $elo_elo=0;

for my $c (keys %call){
    for my $p (keys %{$call{$c}}){
        if ($call{$c}{$p} eq "ano_ano"){ $ano_ano++;}
        if ($call{$c}{$p} eq "ano_tru"){ $ano_tru++;}
        if ($call{$c}{$p} eq "ano_elo"){ $ano_elo++;}
        if ($call{$c}{$p} eq "tru_ano"){ $tru_ano++;}
        if ($call{$c}{$p} eq "tru_tru"){ $tru_tru++;}
        if ($call{$c}{$p} eq "tru_elo"){ $tru_elo++;}
        if ($call{$c}{$p} eq "elo_ano"){ $elo_ano++;}
        if ($call{$c}{$p} eq "elo_tru"){ $elo_tru++;}
        if ($call{$c}{$p} eq "elo_elo"){ $elo_elo++;}
    }
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
print ",,nterminal_support,,\n";
print "predictions,,ann,elo,tru\n";
print ",sum,$n_sum_ano,$n_sum_elo,$n_sum_tru\n";
print "ann,$pred_sum_ano,$ano_ano,$ano_elo,$ano_tru\n";
print "ext,$pred_sum_elo,$elo_ano,$elo_elo,$elo_tru\n";
print "tru,$pred_sum_tru,$tru_ano,$tru_elo,$tru_tru\n";

#summary stats
#True +ve: Supported and predicted ORFs 
#False +ve: Predicted ORFs that disagreed with supported positions.
#False -ve: Supported genes where no ORF was predicted.
#True -ve:  All potential TIS in supported stop2stop regions, that were neither predicted nor supported.

print "\nthere are $count_supported, genes with n-termnial info\n";
print "there are $total_tis, candiate TIS in all supported genes\n";

my $predicted_correct=$ano_ano+$elo_elo+$tru_tru;
my $predicted_incorrect=$ano_elo+$ano_tru+$elo_ano+$elo_tru+$tru_ano+$tru_elo;
my $supported_not_found=$count_supported-$predicted_correct-$predicted_incorrect;
my $correctly_rejected=$total_tis-$predicted_correct-$predicted_incorrect-$supported_not_found;

my $sensitivity=$predicted_correct/($predicted_correct+$supported_not_found);
my $specificity=$correctly_rejected/($predicted_incorrect+$correctly_rejected);
my $precision=$predicted_correct/($predicted_correct+$predicted_incorrect);

print "\ntrue positive: $predicted_correct\n";
print "true negative: $correctly_rejected\n";
print "false positive: $predicted_incorrect\n";
print "false negative: $supported_not_found\n";

print "\nsensitivity: $sensitivity\n";
print "specificity: $specificity\n";
print "positive_predicitive_value: $precision\n";

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub stopToStopFromGTF {

    #regerences to gene and fasta hashes
    my $fasta_ref=$_[0];
    my $gene_info_fwd_ref=$_[1];
    my $gene_info_rev_ref=$_[2];
    my $gene_2_chr_ref=$_[3];

    #to dereferece (makes a new copy og the hash)
    my %gene_info_fwd=%{$gene_info_fwd_ref};
    my %gene_info_rev=%{$gene_info_rev_ref};
    my %gene_2_chr=%{$gene_2_chr_ref};

    my %in_frame_fwd;     #chr,pos,gene,1
    my %in_frame_rev;     #key=chr, pos, gene, value=1;

    for my $gene ( keys %gene_info_fwd ) {
        for my $annotated_start (keys %{$gene_info_fwd{$gene}}){

            #get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match)     
            my $search=1;
            my $count=0;  #set upper limit to 1000 (100)
            my $pos=$annotated_start-1; #-1 because of 0 based fasta offset 

            while ($search){
                $pos=$pos-3;

                #check that the substing is not smaller than the chr!
                if ($pos<0){ last; }

                my $seq=substr($fasta_sequences{$gene_2_chr{$gene}},$pos,3);
                #check for stop codon
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                    $search=0;
                }

                if ($count>=999){ $search=0; }
                $count++;
            }

            #also loop for cds and in frame positions
            my $frameCount=0;
            for ($pos+4 .. $gene_info_fwd{$gene}{$annotated_start}){

                $frameCount++;
                if ($frameCount%3 == 1){
                    $in_frame_fwd{$gene_2_chr{$gene}}{$_}{$gene}=1;
                }
            }
        }
    }

    for my $gene (keys %gene_info_rev){
        for my $annotated_start (keys %{$gene_info_rev{$gene}}){

            #get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match)     
            my $search=1;
            my $count=0;  #set upper limit to 1000     
            my $pos=$annotated_start-1;
            while ($search){ #in the upstream stop codon
                $pos=$pos+3;

                #check that the substing is not bigger or smaller than the chr!
                if ($pos+1>length($fasta_sequences{$gene_2_chr{$gene}})){ last; }

                my $seq=reverse(substr($fasta_sequences{$gene_2_chr{$gene}},$pos-2,3));
                $seq=~tr/ACGTacgt/TGCAtgca/;
                #check for stop codon
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                    $search=0;
                }

               if ($count>=999){ $search=0; }
               $count++;
            }

            #also loop for cds and in frame positions (this is the whole stop2stop region)
            my $frameCount=0;
            #start to end
            for ($gene_info_rev{$gene}{$annotated_start} .. $pos-2) {
                $frameCount++;
                if ($frameCount%3 == 1){
                    $in_frame_rev{$gene_2_chr{$gene}}{$_}{$gene}=1;
                }
            }
        }
    }
    return (\%in_frame_fwd, \%in_frame_rev);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub countTIS {

    #regerences to gene and fasta hashes
    my $fasta_ref=$_[0];
    my $gene_info_fwd_ref=$_[1];
    my $gene_info_rev_ref=$_[2];
    my $gene_2_chr_ref=$_[3];
    my $gene_supp_ref=$_[4];

    #to dereferece (makes a new copy og the hash)
    my %gene_info_fwd=%{$gene_info_fwd_ref};
    my %gene_info_rev=%{$gene_info_rev_ref};
    my %gene_2_chr=%{$gene_2_chr_ref};
    my %supported=%{$gene_supp_ref};

    my $total_potential_TIS=0;

    for my $gene ( keys %gene_info_fwd ) {
        if (exists ($supported{$gene})){  #     #restrict to those genes with N-termnial support
            for my $annotated_start (keys %{$gene_info_fwd{$gene}}){

                #get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match)     
                my $search=1;
                my $count=0;  #set upper limit to 1000 (100)
                my $pos=$annotated_start-1;
                while ($search){  #find upstream inframe stop codon
                    $pos=$pos-3;
                    #check that the substing is not smaller than the chr!
                    if ($pos<0){ last; }
                    my $seq=substr($fasta_sequences{$gene_2_chr{$gene}},$pos,3);
                    #check for stop codon
                    if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                        $search=0;
                    }
                    if ($count>=999){ $search=0; }
                    $count++;
                }

                #also loop for cds and in frame positions
                my $frameCount=0;
                for ($pos+4 .. $gene_info_fwd{$gene}{$annotated_start}){
                   $frameCount++;
                    if ($frameCount%3 == 1){
                        #count the number of inframe near canonical positions here. 
                        #restrict to those genes with N-termnial support
                        my $seq=substr($fasta_sequences{$gene_2_chr{$gene}},$_,3);
                        if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){
                            $total_potential_TIS++;
                        }
                    }
                }
            }
        }
    }

    for my $gene (keys %gene_info_rev){
        if (exists ($supported{$gene})){  #     #restrict to those genes with N-termnial support 
            for my $annotated_start (keys %{$gene_info_rev{$gene}}){

                #get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match)     
                my $search=1; 
                my $count=0;  #set upper limit to 1000     
                my $pos=$annotated_start-1;
                while ($search){ #in the upstream stop codon
                    $pos=$pos+3;
                    #check that the substing is not bigger or smaller than the chr! 
                    if ($pos+1>length($fasta_sequences{$gene_2_chr{$gene}})){ last; }
                    my $seq=reverse(substr($fasta_sequences{$gene_2_chr{$gene}},$pos-2,3));
                    $seq=~tr/ACGTacgt/TGCAtgca/;
                    #check for stop codon
                    if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                        $search=0;
                    }
                    if ($count>=999){ $search=0; }
                    $count++;
                }

                #also loop for cds and in frame positions (this is the whole stop2stop region)
                my $frameCount=0;
                #start to end
                for ($gene_info_rev{$gene}{$annotated_start} .. $pos-2) {
                    $frameCount++;
                    if ($frameCount%3 == 1){
                        my $seq=reverse(substr($fasta_sequences{$gene_2_chr{$gene}},$pos-2,3));
                        $seq=~tr/ACGTacgt/TGCAtgca/;
                        if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){
                            $total_potential_TIS++;
                        }
                    }
                }
            }
        }
    }
    return ($total_potential_TIS);
}
