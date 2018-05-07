#!/usr/bin/perl -w
use strict;

#to do 07/05/2018
#script to take a gft and a set of genes with a set validated start codons. Genes with validated start codons will be shuffled and a new randomly selected in frame potential start codons will be selected (potential start codons are AT*, A*G, *TG).

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $nterm=$ARGV[2];    #bed of N-termini predictions
my $out_gtf=$ARGV[3];  #output gtf (sfter shuffling)

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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#Open n-termini bed
#Extend n-terminal peptides to downstream in frame stop codons
#Store the start and stop positions of n-terminally supported ORFs

my %n_term_fwd; #key=chr, key=1st nt of stop codon, value=1st nt of start codon
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
            my $stop_position=&closest_inframe_downstream_stop_codon_fwd($chr,$start);
            $n_term_fwd{$chr}{$stop_position}=$start; 
        }else{                                                 #reverse cases
            #five_prime=$stop;
            my $stop_position=&closest_inframe_downstream_stop_codon_rev($chr,$stop);
            $n_term_rev{$chr}{$stop_position}=$stop;
        }
    }
}
close (PEP);

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Open genes 
#Find annotated genes that share stop codons with an N-terminal peptide

my $total_TIS_count_in_supported_genes=0;
my $supported_ORF_count=0;
my $supported_matching=0;
my $supported_truncation=0;
my $supported_elongation=0;

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

        if ($class eq "CDS"){
		#should I repat this for exons to keep the GTF consistent?

            if ($dir eq "+"){ 
 
                if (exists ($n_term_fwd{$chr}{$prim3-2})){ #compare the stop codons

                    #get all potential TIS in the stop 2 stop regions
                    #randomly select one of them
                    #print out the start and the new stop

                    #get all potenital TIS
                    my $potential_TIS_ref=&find_potential_TIS_in_stop_2_stop_region_fwd($chr,$prim3-2,$prim5);
                    my @potential_TIS=@{$potential_TIS_ref};
                    #print "fwd,$gene_id,$chr"; for (sort @potential_TIS){ print ",$_"; } print "\n"; 
                    print "$chr\tshuffled\t$class\t".$potential_TIS[rand @potential_TIS]."\t$prim3\t$dir\t$b[6]\t$b[7]\t$b[8]\n";
                    
                #do I also need to print the "exon" line

                }else{
                    print "$_";
                #do I also need to print the "exon" line

                }
            }else{

                if (exists ($n_term_rev{$chr}{$prim5+2})){ #compare the stop codons

                    #get all  potential TIS
                    my $potential_TIS_ref=&find_potential_TIS_in_stop_2_stop_region_rev($chr,$prim5+2,$prim3);
                    my @potential_TIS=@{$potential_TIS_ref};
                    #print "rev,$gene_id,$chr"; for (sort @potential_TIS print ",$_"; } print "\n";
                    print "$chr\tshuffled\t$class\t$prim5\t".$potential_TIS[rand @potential_TIS]."\t$b[6]\t$dir\t$b[7]\t$b[8]\n";
                #do I also need to print the "exon" line
                                }else{
                    print "$_";
                #do I also need to print the "exon" line
                }
            }
        }
    }
}
close(GENES);

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub closest_inframe_downstream_stop_codon_fwd{

    my $chr=$_[0];
    my $start_position=$_[1];

    if (exists ($fasta_sequences{$chr})){

        my $search=1;
        while ($search){ 
            $start_position=$start_position+3;
                    
            #check that the position is not larger than the chromosome
            if ($start_position>((length($fasta_sequences{$chr})-1))){ last; }

            my $seq=substr($fasta_sequences{$chr},($start_position-1),3);
             
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
        }
    }    
    return ($start_position); #the 1st nuclotide of the stop codon
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub closest_inframe_downstream_stop_codon_rev{
    
    my $chr=$_[0];
    my $start_position=$_[1];
    
    if (exists ($fasta_sequences{$chr})){
        
        my $search=1; 
        while ($search){  #find upstream inframe stop codon
            $start_position=$start_position-3;
                    
            #check that the substing is not smaller than the chromosome
            if ($start_position<0){ last; }
            
            my $seq=reverse(substr($fasta_sequences{$chr},$start_position-3,3));
            $seq=~tr/ACGTacgt/TGCAtgca/;
                    
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
        }

    }
    return $start_position; #the 1st nuclotide of the stop codon
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub find_potential_TIS_in_stop_2_stop_region_fwd{

    my $chr=$_[0];
    my $stop_codon=$_[1];
    my $start_codon=$_[2];

    my @potential_TIS;

    if (exists ($fasta_sequences{$chr})){

        #start at stop codon and procees upstream 3nt at a time until we find a stop codon.
        #or have gone more than 999nt past the start codon.
    
        my $search=1;
        my $count=0;
        my $upstream_of_TIS=0;
 
        while ($search){ 
            $stop_codon=$stop_codon-3;

            #check that the substing is not smaller than the chr!
            if ($stop_codon<0){ last; }
 
            my $seq=substr($fasta_sequences{$chr},($stop_codon-1),3);

            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
    
            if ($count>=999){ $search=0; }
          
            if ($upstream_of_TIS) { $count+=3; }

           if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ push (@potential_TIS, $stop_codon); }
#            if ($seq=~/ATG/){ push(@potential_TIS, $stop_codon); }

            if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }

        }
    }
    return(\@potential_TIS);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub find_potential_TIS_in_stop_2_stop_region_rev{

    my $chr=$_[0];
    my $stop_codon=$_[1];
    my $start_codon=$_[2];

    my @potential_TIS;

    if (exists ($fasta_sequences{$chr})){

        #start at stop codon and procees upstream 3nt at a time until we find a stop codon.
        #or have gone more than 999nt past the start codon.
    
        my $search=1;
        my $count=0;  
        my $upstream_of_TIS=0;
 
        while ($search){ 
            $stop_codon=$stop_codon+3;
 
            #check that the substing is not larger than the chr
            if ($stop_codon>=(length($fasta_sequences{$chr})-1)){ last; }
 
            my $seq=reverse(substr($fasta_sequences{$chr},$stop_codon-3,3));
            $seq=~tr/ACGTacgt/TGCAtgca/;
 
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
            
            if ($count>=999){ $search=0; }
            
            if ($upstream_of_TIS) { $count+=3; }

            if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ push(@potential_TIS, $stop_codon); }
#            if ($seq=~/ATG/){ push(@potential_TIS, $stop_codon); }

            if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }

        }
    }
    return(\@potential_TIS);
    #this should be a reference to the array
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
