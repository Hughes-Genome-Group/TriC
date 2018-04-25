#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#description:
#this script performs in silico digestion of a fastq file and outputs the digested reads >20 bp in fastq format
#the number identifying the pair-mate (1 or 2) and restriction fragment (0 to n) is reported after the id line in the fastq file
#the restriction enzyme sequence needs to be specified by the user

#bugs:
#if the read contains more than 3 restriction enzyme cut sites in a row at the start or end of the read, the read won't be reported (but this will hardly ever/never happen)

#run:
#run script with full path to file to be digested, eg: perl fastq_digester.pl /t1-data1/WTSA_Dev/oudelaar/triC/pilot/Rep1_combined_reads.fastq

#to be specified by user
my $RE_seq = "CATG";                                                        #restriction enzyme sequence, in this case NlaIII

#file names
my $fastq_file = $ARGV[0];
my $path = "undef";
my $file = "undef";
if ($fastq_file =~ /(.*\/)(.*)\.(fastq|fq)/) {
    $path = $1;
    $file = $2;
    }
my $fastq_out = "$path$file\_dig.fastq";

#initiate variables
my %fq_hash;
my @line_labels = ("id", "seq", "plus", "qscore");
my $line_counter = 0;
my $seq_counter = 0;
my $min_length = 20;        

#open file handles
open (FQ_FH, $fastq_file) or die "can't open fastq file $fastq_file";
open (OUT, ">$fastq_out") or die "can't open fastq output file $fastq_out";

#digest fastq files
while ($fq_hash{$line_labels[$seq_counter]} = <FQ_FH>) {                    #read fastq file, 4 lines at the time, store in hash 
    chomp $fq_hash{$line_labels[$seq_counter]};
    $line_counter++;
    $seq_counter++;
    if ($seq_counter == 4) {
        my $id_begin = "undef";
        my $id_end = "undef";
        my $PE_read = "undef";
        if ($fq_hash{"id"} =~ /(.*) (\d):(.*:.*:.*)/ ) {
            $id_begin = $1;
            $id_end = " ".$2.":".$3;
            $PE_read = $2;
        }
        my @rf_array = split(/$RE_seq/, $fq_hash{"seq"});                   #split sequence in separate restriction fragments
        my $q_tracker = 0;                                                  #sum length of fragments while going through the loop to output corresponding q-score
        for (my $i = 0; $i < $#rf_array + 1; $i++) {
            $q_tracker = length($rf_array[$i]) + $q_tracker + 4;
            if ($#rf_array == 0 and $fq_hash{"seq"} =~ /^$RE_seq/) {        
                if ((length($rf_array[$i]) + 4) == length($fq_hash{"qscore"})) {        #restriction site at start of sequence: 1 site
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$RE_seq$rf_array[$i]\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                elsif ((length($rf_array[$i]) + 8) == length($fq_hash{"qscore"})) {     #restriction site at start of sequence: 2 sites
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$RE_seq$RE_seq$rf_array[$i]\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                elsif ((length($rf_array[$i]) + 12) == length($fq_hash{"qscore"})) {    #restriction site at start of sequence: 3 sites
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$RE_seq$RE_seq$RE_seq$rf_array[$i]\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                else {
                    print "$fq_hash{\"seq\"}\n";
                }
            }
            elsif ($#rf_array == 0 and $fq_hash{"seq"} =~ /$RE_seq$/) {     
                if ((length($rf_array[$i]) + 4) == length($fq_hash{"qscore"})) {        #restriction site at end of sequence: 1 site
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$rf_array[$i]$RE_seq\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                elsif ((length($rf_array[$i]) + 8) == length($fq_hash{"qscore"})) {     #restriction site at end of sequence: 2 sites
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$rf_array[$i]$RE_seq$RE_seq\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                elsif ((length($rf_array[$i]) + 12) == length($fq_hash{"qscore"})) {    #restriction site at end of sequence: 3 sites
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$rf_array[$i]$RE_seq$RE_seq$RE_seq\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
                }
                else {
                    print "$fq_hash{\"seq\"}\n";
                }
            }
            elsif ($#rf_array == 0) {                                       #no restriction site
                print OUT "$id_begin:RF$PE_read:$i $id_end\n$rf_array[$i]\n$fq_hash{\"plus\"}\n$fq_hash{\"qscore\"}\n";
            }
            elsif (length($rf_array[$i]) >= $min_length - 8) { 
                if ($i == 0) {                                              #first fragment
                    my $q = substr($fq_hash{"qscore"}, 0, length($rf_array[$i]) + 4);
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$rf_array[$i]$RE_seq\n$fq_hash{\"plus\"}\n$q\n";
                }
                elsif ($i == $#rf_array) {                                  #last fragment
                    my $q = substr($fq_hash{"qscore"}, -length($rf_array[$i]) - 4, $q_tracker);          
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$RE_seq$rf_array[$i]\n$fq_hash{\"plus\"}\n$q\n";
                }
                else {                                                      #middle fragments
                    my $q = substr($fq_hash{"qscore"}, $q_tracker - length($rf_array[$i]) - 8, length($rf_array[$i]) + 8);
                    print OUT "$id_begin:RF$PE_read:$i $id_end\n$RE_seq$rf_array[$i]$RE_seq\n$fq_hash{\"plus\"}\n$q\n";
                }   
            }
        }
        $seq_counter = 0;
    }
}


