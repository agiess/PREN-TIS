#!/usr/bin/perl -w
use strict;

#to do 22/01/2018
#script to produce the matrix for pca, start codon classifiaction
#output is split into lengths.
#for each gene look upstream of the start codon to find the next in frame stop codon
#splitting codon counts per codon

#with fpkm 
#whole genome

#input files:
#gtf exons, start, stop codons
#sam = aligned 
my $gtf=$ARGV[0];
my $sam=$ARGV[1];
my $fasta=$ARGV[2];
my $out_file=$ARGV[3];
my $min_length=$ARGV[4];
my $max_length=$ARGV[5];

my $WINDOW_START=20; #20 nt upstream
my $WINDOW_END=20;   #20 nt downstream

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
my %start_codon; #key=gene_id, value=pos
my %gene_info_fwd; #key=gene_id, key2=start, value=stop
my %gene_info_rev; #key=gene_id, key2=start, value=stop
my %gene_lengths; #for FPKM
my %gene_2_chr; #key = gene_id; value = chr

my %ann_start_codon_fwd; #key=chr, key2=pos, value=1   #used for flags only
my %ann_start_codon_rev; #key=chr, key2=pos, value=1   
my %multi_exon;

open(GENES,$gtf) || die "can't open $gtf";        #gtf is 1 based
while (<GENES>){
    unless(/^#/){
        my @b=split("\t");
        my $class=$b[2];
        my $chr=$b[0];
        my $start=$b[3];
        my $end=$b[4];
        my $dir=$b[6];
        my $gene_id="unknown";
        ($gene_id) = $b[8] =~ /gene_id\s"([^\"]+)";/;

        if ($class eq "CDS"){ #reworked to account for gtf's without start/stop codons

            if ($dir eq "+"){ #fwd cases. Use start positions as 5'

     #           #fwd genes take start stop from first exons.

                if (exists( $multi_exon{$gene_id}) ){
                    $gene_lengths{$gene_id}+=($end+1)-$start;  #add on this start to end
                    my $previous_start=$start_codon{$gene_id};
                    $gene_info_fwd{$gene_id}{$previous_start}=($end-2); #stops  
                }else{
                    $start_codon{$gene_id}=$start; #starts
                    $ann_start_codon_fwd{$chr}{$start}=1;
                    $gene_info_fwd{$gene_id}{$start}=($end-2); #stops      
                    $gene_lengths{$gene_id}=($end+1)-$start;  #start to end
                    $gene_2_chr{$gene_id}=$chr;
                    $multi_exon{$gene_id}=1
                }

           }else{  #reverse cases. Use stop positions as 5'

    #           #rev genes take start and stop from first exon

                if (exists ( $multi_exon{$gene_id}) ){
                    $gene_lengths{$gene_id}+=($end+1)-$start; #add on this end to start
                    my $previous_start=$start_codon{$gene_id};
                    $gene_info_rev{$gene_id}{$previous_start}=($start+2);
                }else{
                    $start_codon{$gene_id}=$end;
                    $ann_start_codon_rev{$chr}{$end}=1;
                    $gene_info_rev{$gene_id}{$end}=($start+2);
                    $gene_lengths{$gene_id}=($end+1)-$start; #end to start
                    $gene_2_chr{$gene_id}=$chr;
                    $multi_exon{$gene_id}=1;
                }
            }
        }
    }
}
close(GENES);

print "5' gtf parsed\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#open fasta for codon sequeces
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

print "fasta parsed\n";

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#assign regions
my %stop_2_stop_fwd;
my %stop_2_stop_rev;

for my $gene (keys %gene_info_fwd){
    for my $annotated_start (keys %{$gene_info_fwd{$gene}}){ 
        my $stop_codon=$gene_info_fwd{$gene}{$annotated_start};

		#get first nucleotide of start codon then procees upstream 3nt at a time until we find a stop codon (pattern match) 	
        my $search=1; 
        my $count=0;  #set upper limit to 999
        my $pos=$stop_codon; #-4 to skip the stop codon

        while ($search){
            $pos=$pos-3;   

            #check that the substing is not smaller than the chr!
            if ($pos-3 < 0){ $search=0; }
                         
            $stop_2_stop_fwd{ $gene_2_chr{$gene} } {$pos}=1;
                   
            #check for stop codon
            my $seq=substr($fasta_sequences{$gene_2_chr{$gene}},$pos-1,3);  
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }

            if ($count>=333){ $search=0; }
          
            #extend up to 999nt past the start codon
            if ($pos < $annotated_start){
                $count++;
            }
        }
    }
}

for my $gene (keys %gene_info_rev){
    for my $annotated_start (keys %{$gene_info_rev{$gene}}){
        my $stop_codon=$gene_info_rev{$gene}{$annotated_start};

        #get first nucleotide of start codon then proceed upstream 3nt at a time until we find a stop codon (pattern match)     
        my $search=1;
        my $count=0;  #set upper limit to 999     
        my $pos=$stop_codon;

        while ($search){
            $pos=$pos+3;

            #check that the substing is not bigger or smaller than the chr!
            if ($pos+3>length($fasta_sequences{$gene_2_chr{$gene}})){ $search=0 }

            $stop_2_stop_rev{ $gene_2_chr{$gene} } {$pos}=1;

            #check for stop codon
            my $seq=reverse(substr($fasta_sequences{$gene_2_chr{$gene}},$pos-3,3)); #(-2 to leftmost codon position, -1 for fastq offset)
            $seq=~tr/ACGTacgt/TGCAtgca/;
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
            if ($count>=333){ $search=0; }

            #exetend up to 999nt past the start codon
 
           if ($pos > $annotated_start){
                $count++;
            }
        }
    } 
}

###################################################################################################
#open sam file, store riboseq counts per position
my %counts_fwd; #key = chr, key2 = position, key3 = lenghts, value = counts.
my %counts_rev; #key = chr, key2 = position, key3 = lenghts, value = counts.
my $total_counts=0; #for FPKM 

#If minimum and maxmium read lengths were not provided, calculate them from the sam file 
my $LONGEST_FRAGMENT=1;             #Default_values
my $SHORTEST_FRAGMENT=1000000;      #Default_values

my $find_short=0;
my $find_long=0;

unless ($min_length){
   $find_short=1;
}else{
   $SHORTEST_FRAGMENT=$min_length;
}

unless ($max_length){
   $find_long=1;
}else{
   $LONGEST_FRAGMENT=$max_length;
}

open (SAM, $sam) || die "can't open $sam";      #sam is 1 based
while (<SAM>){
    #skip headers
    unless(/^#/){
        my @b=split("\t");
        my $leftMost=$b[3]; #leftmost position of match 5' for fwd, 3' for rev
        my $flag=$b[1];
        my $chr=$b[2];
        my $mapq=$b[4];
        my $cig=$b[5];
        my $seq=$b[9];
        my $fivePrime;                     
        my $length=length($seq);                        
        
        if ($mapq >= 10){ #mapping uniqness filter
            unless ($flag & 0x4){   #if aligned

                $total_counts++;

                if ($flag & 0x10){  #if rev calculate 5' from the rightmost position

                    #parse cigar for indels and adjust the length of the alignment             
                    my $length5=length($seq);
                    while ($cig =~/(\d+)I/g){   #add to length for insertions
                        $length5+=$1;    
                    }
                    while ($cig =~/(\d+)D/g){   #substact from length for deletions
                        $length5-=$1;
                    }

                    $fivePrime=$leftMost+$length5-1;   #sam is 1 based
                    $counts_rev{$chr}{$fivePrime}{$length}++;

                }else{ #if fwd this is easy
                    $fivePrime=$leftMost;
                    $counts_fwd{$chr}{$fivePrime}{$length}++;
                }

				#find the fragment length range for the output matrix 
                if ($find_short){
                    if ($length<$SHORTEST_FRAGMENT){
				        $SHORTEST_FRAGMENT=$length;
                    }
                }
 
                if ($find_long){
                    if ($length>$LONGEST_FRAGMENT){
                        $LONGEST_FRAGMENT=$length;
                    }
				}
            }
        }else{
            unless ($flag & 0x4){   #if aligned
                $total_counts++;
            }   
        }
    }                    
}
close (SAM);

print "sam parsed\n";

####################################################################################################
#summarise and output matrix
open (OUT, ">$out_file") || die;

#header
#print OUT "#id,codon,dir,canonical_candidate_sites_ORF_region,annotated_start_site,near_cognate_codon,reads_at_pos,window_reads_downstream,window_reads_upstream,proportion_of_reads_at_position,proportion_of_reads_upstream,proportion_of_reads_downstream,ORF_fpkm";
#print OUT "#id,codon,dir,distance_to_region_start,annotated_start_site,near_cognate_codon,reads_at_pos,window_reads_downstream,window_reads_upstream,proportion_of_reads_at_position,proportion_of_reads_upstream,proportion_of_reads_downstream,ORF_fpkm";
print OUT "#id,codon,dir,upstream_ATG_count,upstream_AAG_count,upstream_ACG_count,upstream_AGG_count,upstream_CTG_count,upstream_GTG_count,upstream_TTG_count,upstream_ATA_count,upstream_ATC_count,upstream_ATT_count,annotated_start_site,near_cognate_codon,reads_at_pos,window_reads_downstream,window_reads_upstream,proportion_of_reads_at_position,proportion_of_reads_upstream,proportion_of_reads_downstream,ORF_fpkm";

for my $hpos (-$WINDOW_START .. $WINDOW_END){
    print OUT ",seq_".$hpos;    
}

for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
	for my $pos(-$WINDOW_START .. $WINDOW_END){
         print OUT ",".$pos."_".$len;
	}
}
 
print OUT "\n";

#Parse once for reverse stand
for my $chr (sort keys %fasta_sequences){

    my $end=(length($fasta_sequences{$chr}));  #rightmost  #zero or 1 based?
    my $start=0;                               #leftmost

    #windows
    for my $codon_pos ($start+$WINDOW_END .. $end-$WINDOW_START){
        my $annotated_start_site="FALSE";
        my $near_cognate="FALSE";
        my $reads_at_pos=0;
        my $window_reads_downstream=0;
        my $window_reads_upstream=0;
        my $win_sum=0;

        #take inframe codons in stop to stop regions
        if (exists ($stop_2_stop_rev{$chr}{$codon_pos})){ 

            my $codon=reverse(substr($fasta_sequences{$chr},$codon_pos-3,3));    
            $codon=~tr/ACGTacgt/TGCAtgca/;

            if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){ #take only cognate and near cognate codons
                $near_cognate="TRUE";
                if (exists ($ann_start_codon_rev{$chr}{$codon_pos})){$annotated_start_site="TRUE";}
 
                #parse all lenghts for count summaries (summmed over all fragment lengths)
                for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                    for my $win_pos (($codon_pos-$WINDOW_END) .. ($codon_pos+$WINDOW_START)){
                        # up and downstream should be opposite for rev reads
                        if (exists ($counts_rev{$chr}{$win_pos}{$len})){
                            $win_sum+=$counts_rev{$chr}{$win_pos}{$len};
                            if ($win_pos==$codon_pos){ $reads_at_pos+=$counts_rev{$chr}{$win_pos}{$len}; }
                            if ($win_pos<$codon_pos){ $window_reads_downstream+=$counts_rev{$chr}{$win_pos}{$len}; }
                            if ($win_pos>$codon_pos){ $window_reads_upstream+=$counts_rev{$chr}{$win_pos}{$len}; }
                        }
                    }
                }

                my $proportion_at_position=eval {$reads_at_pos/$win_sum} || 0;
                my $proportion_upstream=eval {$window_reads_upstream/$win_sum} || 0;
                my $proportion_downstream=eval {$window_reads_downstream/$win_sum} || 0;

                my $id=$chr."_".$codon_pos."_rev";
                my $window_seq=reverse(substr($fasta_sequences{$chr}, ($codon_pos-($WINDOW_START+1)), ($WINDOW_START+$WINDOW_END+1)));  
                $window_seq=~tr/ACGTacgt/TGCAtgca/;
 
                #fkpm + codon rank here
                #my ($ORF_FPKM, $codon_rank)=&stop2stop_rev($chr,$codon_pos);
                #my ($ORF_FPKM, $distance_to_start_of_stop2_stop_region)=&stop2stop_rev($chr,$codon_pos);

                my ($ORF_FPKM, $upstream_ATG, $upstream_AAG, $upstream_ACG, $upstream_AGG, $upstream_CTG,$upstream_GTG, $upstream_TTG, $upstream_ATA, $upstream_ATC, $upstream_ATT)=&stop2stop_rev($chr,$codon_pos);
              
                #print OUT "$id,$codon,rev,$codon_rank,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";
                #print OUT "$id,$codon,rev,$distance_to_start_of_stop2_stop_region,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";
              
                print OUT "$id,$codon,rev,$upstream_ATG,$upstream_AAG,$upstream_ACG,$upstream_AGG,$upstream_CTG,$upstream_GTG,$upstream_TTG,$upstream_ATA,$upstream_ATC,$upstream_ATT,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";  

                #Loop through sequence here
                my @seq_out=split("",$window_seq); #window seq is already reversed
                for my $nuc (@seq_out){
                    print OUT ",$nuc";
                }

                #LOOP THROUGH LENGTHS HERE
                for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                    my $win_pos=$codon_pos+$WINDOW_START;
                    while ($win_pos >= $codon_pos-$WINDOW_END){
                        my $signal=0;
                        #riboseq signal       #%counts; #key = chr, key2 = position, value = counts.
                        if (exists ($counts_rev{$chr}{$win_pos}{$len})){
                             $signal=$counts_rev{$chr}{$win_pos}{$len};
                        }
                   
                        my $fraction_of_window_signal=eval {$signal/$win_sum} || 0;
                        print OUT ",$fraction_of_window_signal";
                        $win_pos--;
                    } 
                }
                print OUT "\n";
            }
        }
    }
}

print "rev genes processed\n";

#loop though chromosomes a second time for forwards genes;
for my $chr (sort keys %fasta_sequences){

    my $end=(length($fasta_sequences{$chr}));  #rightmost
    my $start=0;                               #leftmost

    #windows   
    for my $codon_pos ($start+$WINDOW_START .. $end-$WINDOW_END){
        my $annotated_start_site="FALSE";
        my $near_cognate="FALSE";
        my $reads_at_pos=0;
        my $window_reads_downstream=0;
        my $window_reads_upstream=0;
        my $win_sum=0;

        #take inframe codons in stop to stop regions
        if (exists ($stop_2_stop_fwd{$chr}{$codon_pos})){ 
  
            my $codon=substr($fasta_sequences{$chr},$codon_pos-1,3);
            if ($codon =~/[ACGT]TG/ || $codon =~/A[ACGT]G/ || $codon =~/AT[ACGT]/){ #take only congate and near cognate codons

                $near_cognate="TRUE";
                if (exists ($ann_start_codon_fwd{$chr}{$codon_pos})){$annotated_start_site="TRUE";}
        
                #parse all lenghts for count summaries (summmed over all fragment lengths)
                for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                    for my $win_pos (($codon_pos-$WINDOW_START) .. ($codon_pos+$WINDOW_END)){
                        if (exists ($counts_fwd{$chr}{$win_pos}{$len})){
                            $win_sum+=$counts_fwd{$chr}{$win_pos}{$len};
                            if ($win_pos==$codon_pos){ $reads_at_pos+=$counts_fwd{$chr}{$win_pos}{$len}; }
                            if ($win_pos>$codon_pos){ $window_reads_downstream+=$counts_fwd{$chr}{$win_pos}{$len}; }
                            if ($win_pos<$codon_pos){ $window_reads_upstream+=$counts_fwd{$chr}{$win_pos}{$len}; }
                        }
                    }
                }
                
                #my $window_up_down_ratio=eval {$window_reads_downstream/$window_reads_upstream} || 0;
                my $proportion_at_position=eval {$reads_at_pos/$win_sum} || 0;
                my $proportion_upstream=eval {$window_reads_upstream/$win_sum} || 0;
                my $proportion_downstream=eval {$window_reads_downstream/$win_sum} || 0;

                my $id=$chr."_".$codon_pos."_fwd";
                my $window_seq=substr($fasta_sequences{$chr}, ($codon_pos-($WINDOW_START+1)), ($WINDOW_START+$WINDOW_END+1));
 
                #fkpm + codon rank here
                #my ($ORF_FPKM,$codon_rank)=&stop2stop_fwd($chr,$codon_pos);
                #my ($ORF_FPKM,$distance_to_start_of_stop2_stop_region)=&stop2stop_fwd($chr,$codon_pos);
                my ($ORF_FPKM, $upstream_ATG, $upstream_AAG, $upstream_ACG, $upstream_AGG, $upstream_CTG, $upstream_GTG, $upstream_TTG, $upstream_ATA, $upstream_ATC, $upstream_ATT)=&stop2stop_fwd($chr,$codon_pos);

                #print OUT "$id,$codon,fwd,$codon_rank,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";
                #print OUT "$id,$codon,fwd,$distance_to_start_of_stop2_stop_region,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";
      
                print OUT "$id,$codon,fwd,$upstream_ATG,$upstream_AAG,$upstream_ACG,$upstream_AGG,$upstream_CTG,$upstream_GTG,$upstream_TTG,$upstream_ATA,$upstream_ATC,$upstream_ATT,$annotated_start_site,$near_cognate,$reads_at_pos,$window_reads_downstream,$window_reads_upstream,$proportion_at_position,$proportion_downstream,$proportion_upstream,$ORF_FPKM";

                #Loop through sequence here
                #output sequence
                my @seq_out=split("",$window_seq);
                for my $nuc (@seq_out){
                    print OUT ",$nuc";
                }

                #LOOP THROUGH LENGTHS HERE
                for my $len ($SHORTEST_FRAGMENT .. $LONGEST_FRAGMENT){
                    #second loop for output
                    for my $win_pos (($codon_pos-$WINDOW_START) .. ($codon_pos+$WINDOW_END)){                  
                        my $signal=0;
                        #riboseq signal
                        #%counts; #key = chr, key2 = position, value = counts.
                        if (exists ($counts_fwd{$chr}{$win_pos}{$len})){ 
                            $signal+=$counts_fwd{$chr}{$win_pos}{$len};
                        }

                        my $fraction_of_window_signal=eval {$signal/$win_sum} || 0;            
                            print OUT ",$fraction_of_window_signal";
                    }
                }
                print OUT "\n";
            }
        }
    }
}

close (OUT);

print "fwd gene processed\nshortest: $SHORTEST_FRAGMENT\nlongest: $LONGEST_FRAGMENT\n";

exit;

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub stop2stop_fwd{
    #for a given position (chr + pos)
    #find the nearest downstream in frame stop codon, and upstream in frame stop codon (to a limit of 999nt, without going out of chr limits)
    #calculate fkpm
    #calculate the number of potential start sites upstream of this one + 1

    my $CHR=$_[0];
    my $current_codon_pos=$_[1];
    my $upstream_codons=0;
    my $distance_to_upstream_stop_codon=0;

    my $upstream_ATG=0;
    my $upstream_AAG=0;
    my $upstream_ACG=0;
    my $upstream_AGG=0;
    my $upstream_CTG=0;
    my $upstream_GTG=0;
    my $upstream_TTG=0;
    my $upstream_ATA=0;
    my $upstream_ATC=0;
    my $upstream_ATT=0;

    ###
    #find upstream stop codon
    ###

    #keep track of the potential start sites     
    my $search=1;
    my $count=0;  #set upper limit to 999
    my $pos=$current_codon_pos-1;
    while ($search){
        $pos=$pos-3;

        #check that the substing is not smaller than the chr!
        if ($pos<0){ last; }

        my $seq=substr($fasta_sequences{$CHR},$pos,3);
        if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
            $search=0;
        }

        if ($seq=~ /[ACGT]TG/ || $seq =~ /A[ACGT]G/ || $seq =~ /AT[ACGT]/ ){
            $upstream_codons++;
        }

        if ($seq=~ /^ATG/){ $upstream_ATG++; }
        if ($seq=~ /^AAG/){ $upstream_AAG++; }
        if ($seq=~ /^ACG/){ $upstream_ACG++; }
        if ($seq=~ /^AGG/){ $upstream_AGG++; }
        if ($seq=~ /^CTG/){ $upstream_CTG++; }
        if ($seq=~ /^GTG/){ $upstream_GTG++; }
        if ($seq=~ /^TTG/){ $upstream_TTG++; }
        if ($seq=~ /^ATA/){ $upstream_ATA++; }
        if ($seq=~ /^ATC/){ $upstream_ATC++; }
        if ($seq=~ /^ATT/){ $upstream_ATT++; }

        if ($count>=333){ $search=0; } #limit to 999nt
        $count++
    }

    $distance_to_upstream_stop_codon=$count*3;

    ###
    #find downstream stop codon
    ###

    $search=1;
    my $limit=length($fasta_sequences{$CHR});  #set upper limit to the end fo the chromosome
    $pos=$current_codon_pos-1;
    my $offset=$pos%3;
    my $ORF_length=0;
    my $read_sum=0;
    while ($search){
        $pos++;
        $ORF_length++;

        #sum reads for fpkm
        if (exists ($counts_fwd{$CHR}{$pos})){
            for my $length (keys %{ $counts_fwd{$CHR}{$pos} } ){
                $read_sum+=$counts_fwd{$CHR}{$pos}{$length};
            }
        }

        if ($pos%3 == $offset){

            #check that the substing is not smaller than the chr!
            if ($pos+3>$limit){ $search=0; }

            my $seq=substr($fasta_sequences{$CHR},$pos,3);
            #check for stop codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
               $search=0;
            }
        }
    }

    my $ORF_FPKM=eval { (1000000000*$read_sum)/($total_counts*$ORF_length) } || 0;
    #return ($ORF_FPKM, $upstream_codons);
    #return ($ORF_FPKM, $distance_to_upstream_stop_codon);
    return($ORF_FPKM, $upstream_ATG, $upstream_AAG, $upstream_ACG, $upstream_AGG, $upstream_CTG,$upstream_GTG, $upstream_TTG, $upstream_ATA, $upstream_ATC, $upstream_ATT);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
sub stop2stop_rev{

    my $CHR=$_[0];
    my $current_codon_pos=$_[1];

    my $upstream_codons=0;
    my $ORF_stop=0;
    my $distance_to_upstream_stop_codon=0;
 
    my $upstream_ATG=0;
    my $upstream_AAG=0;
    my $upstream_ACG=0;
    my $upstream_AGG=0;
    my $upstream_CTG=0;
    my $upstream_GTG=0;
    my $upstream_TTG=0;
    my $upstream_ATA=0;
    my $upstream_ATC=0;
    my $upstream_ATT=0;

    ###
    #find upstream stop codon
    ###

    my $search=1;
    my $count=0;  #set upper limit to 999 
    my $pos=$current_codon_pos-1;
    while ($search){
        $pos=$pos+3;

        #check that the substing is not bigger or smaller than the chr!
        if ($pos>length($fasta_sequences{$CHR})){ last; }

        my $seq=reverse(substr($fasta_sequences{$CHR},$pos-2,3));
        $seq=~tr/ACGTacgt/TGCAtgca/;
        #check for start codon
        if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
            $search=0;
        }

        if ($seq=~ /[AGCT]TG/ || $seq =~ /A[ACTG]G/ || $seq =~ /AT[ACGT]/ ){
            $upstream_codons++;
        }

        if ($seq=~ /^ATG/){ $upstream_ATG++; }
        if ($seq=~ /^AAG/){ $upstream_AAG++; }
        if ($seq=~ /^ACG/){ $upstream_ACG++; }
        if ($seq=~ /^AGG/){ $upstream_AGG++; }
        if ($seq=~ /^CTG/){ $upstream_CTG++; }
        if ($seq=~ /^GTG/){ $upstream_GTG++; }
        if ($seq=~ /^TTG/){ $upstream_TTG++; }
        if ($seq=~ /^ATA/){ $upstream_ATA++; }
        if ($seq=~ /^ATC/){ $upstream_ATC++; }
        if ($seq=~ /^ATT/){ $upstream_ATT++; }

        if ($count>=333){ $search=0; } #limit to 999nt
        $count++;        
    }

    $distance_to_upstream_stop_codon=$count*3;

    ###
    #find downstream stop codon
    ###

    #get first nucleotide of start codon then proceed upstream 3nt at a time until we find a stop codon (pattern match)     
    $search=1;
    my $limit=0;  #set upper limit to start of chromosome
    $pos=$current_codon_pos-1;
    my $offset=$pos%3;
    my $ORF_length=0;
    my $read_sum=0;
    while ($search){
        $pos--;
        $ORF_length++;

        #sum reads for fpkm
        if (exists ($counts_rev{$CHR}{$pos+2})){
            for my $length (keys %{ $counts_rev{$CHR}{$pos+2} } ){
                $read_sum+=$counts_rev{$CHR}{$pos+2}{$length};
            }
        }

        if ($pos%3 == $offset){   

            #check that the substring is not bigger or smaller than the chr!
            if ($pos-3<=$limit){ $search=0;}

            my $seq=reverse(substr($fasta_sequences{$CHR},$pos-2,3));
            $seq=~tr/ACGTacgt/TGCAtgca/;
            #check for start codon
            if ($seq=~/TAG/ || $seq=~/TAA/ || $seq=~/TGA/ ){
                $search=0;
            }
        }
    }

    my $ORF_FPKM=eval { (1000000000*$read_sum)/($total_counts*$ORF_length) } || 0;  
    #return ($ORF_FPKM, $upstream_codons);
    #return ($ORF_FPKM, $distance_to_upstream_stop_codon);
    return($ORF_FPKM, $upstream_ATG, $upstream_AAG, $upstream_ACG, $upstream_AGG, $upstream_CTG, $upstream_GTG, $upstream_TTG, $upstream_ATA, $upstream_ATC, $upstream_ATT);
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
