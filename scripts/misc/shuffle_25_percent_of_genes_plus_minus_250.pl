#!/usr/bin/perl -w
use strict;

#to do 07/05/2018
#script to take a gft and a set of genes with a set validated start codons. Genes with validated start codons will be shuffled and a new randomly selected in frame potential start codons will be selected (potential start codons are AT*, A*G, *TG).

my $gtf=$ARGV[0];
my $fasta=$ARGV[1]; 
my $out_gtf=$ARGV[2];   #output gtf (sfter shuffling)
my $out_psudo=$ARGV[3]; #output the original gene coord to a psudo N-termini file (to compaer against)

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
#Open genes 

open (OUT1, ">$out_gtf") || die "can't open $out_gtf"; 
open (OUT2, ">$out_psudo") || die "can't open $out_gtf";

open(GENES,$gtf) || die "can't open $gtf";        #gtf is 1 based
while (<GENES>){

    chomp;

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
   
           my $random_number = int(1 + rand(5 - 1));  #generate a random number between 1 and 4
 
           if ($dir eq "+"){ 
 
                if ($random_number == 4){ #reassign #and print to psudo N-termini file
 
                    #print in bed format
                    print OUT2 "$chr\t".($b[3]-1)."\t$b[4]\t$gene_id\t1\t+\n";

                    #get all potential TIS in the stop 2 stop regions
                    #randomly select one of them
                    #print out the start and the new stop
                    #print the original ORF to the psudo N-terminal file

                    #get all potenital TIS
                    my $potential_TIS_ref=&find_potential_TIS_in_stop_2_stop_region_fwd($chr,$prim3-2,$prim5);
                    my @potential_TIS=@{$potential_TIS_ref};
                    #print "fwd,$gene_id,$chr"; for (sort @potential_TIS){ print ",$_"; } print "\n";
                    my $selected_TIS=$potential_TIS[rand @potential_TIS];
                    print OUT1 "$chr\tshuffled\t$class\t$selected_TIS\t$prim3\t$dir\t$b[6]\t$b[7]\t$b[8]\n";
                    print OUT1 "$chr\tshuffled\texon\t$selected_TIS\t$prim3\t$dir\t$b[6]\t$b[7]\t$b[8]\n";
                }else{
                    print OUT1 "$_\n";
                    print OUT1 "$b[0]\t$b[1]\texon\t$b[3]\t$b[4]\t$b[5]\t$b[6]\t$b[7]\t$b[8]\n";
                }
            }else{

                if ($random_number == 4){ #reassign #and print to psudo N-termini file

                    print OUT2 "$chr\t".($b[3]-1)."\t$b[4]\t$gene_id\t1\t-\n";

                    #get all  potential TIS
                    my $potential_TIS_ref=&find_potential_TIS_in_stop_2_stop_region_rev($chr,$prim5+2,$prim3);
                    my @potential_TIS=@{$potential_TIS_ref};
                    #print "rev,$gene_id,$chr"; for (sort @potential_TIS print ",$_"; } print "\n";
                    my $selected_TIS=$potential_TIS[rand @potential_TIS];                   
                    print OUT1 "$chr\tshuffled\t$class\t$prim5\t$selected_TIS\t$b[6]\t$dir\t$b[7]\t$b[8]\n";
                    print OUT1 "$chr\tshuffled\texon\t$prim5\t$selected_TIS\t$b[6]\t$dir\t$b[7]\t$b[8]\n";
                }else{
                    print OUT1 "$_\n";
                    print OUT1 "$b[0]\t$b[1]\texon\t$b[3]\t$b[4]\t$b[5]\t$b[6]\t$b[7]\t$b[8]\n";
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

            if ( ($stop_codon > ($start_codon-250)) && ($stop_codon < ($start_codon+250)) ){
 
                my $seq=substr($fasta_sequences{$chr},($stop_codon-1),3);

                #check for stop codon
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
    
                if ($count>=999){ $search=0; }
          
                if ($upstream_of_TIS) { $count+=3; }

#                if ($seq=~/ATG/ || $seq=~/GTG/ || $seq=~/TTG/ ){ push (@potential_TIS, $stop_codon); }
                if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ push(@potential_TIS, $stop_codon); }
#                if ($seq=~/ATG/){ push(@potential_TIS, $stop_codon); }

                if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }
            }
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

            if ( ($stop_codon > ($start_codon-250)) && ($stop_codon < ($start_codon+250)) ){
 
                my $seq=reverse(substr($fasta_sequences{$chr},$stop_codon-3,3));
                $seq=~tr/ACGTacgt/TGCAtgca/;
 
                #check for stop codon
                if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){ $search=0; }
            
                if ($count>=999){ $search=0; }
            
                if ($upstream_of_TIS) { $count+=3; }

#                if ($seq=~/ATG/ || $seq=~/GTG/ || $seq=~/TTG/ ){ push(@potential_TIS, $stop_codon); }
                if ($seq=~/\wTG/ || $seq=~/A\wG/ || $seq=~/AT\w/ ){ push(@potential_TIS, $stop_codon); }
#                if ($seq=~/ATG/){ push(@potential_TIS, $stop_codon); }

                if ($stop_codon == $start_codon) { $upstream_of_TIS=1; }
            }
        }
    }
    return(\@potential_TIS);
    #this should be a reference to the array
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
