#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;

#This script converts a file containing the detected interacting restriction fragments per Tri-C read (as outputted by tric_reporters.pl)
#for a viewpoint of interest into a file that reports the interaction frequency between those fragments in the format:
#1:1-100    1:100-200   10

#Run command:
#perl tric_frequencies.pl myfragfile.txt

###############################################################################################################################################################################

#Open file and generate output file (in same directory)
my $full_filename = $ARGV[0];
unless ($full_filename =~ /(.*)(_fragments)\.(.*)/) {die "filename does not match format"};
my $file_name = $1;
my $file_path = "";
if ($file_name =~ /(.*\/)(\V++)/) {$file_path = $1; $file_name = $2};
my $output_filename = "$file_path/$file_name\_freq.txt";

open(FH, $full_filename) or die "cannot open $full_filename";
open(FHOUT, ">$output_filename") or die "cannot open $output_filename";

#Store frequency of interactions between fragments in hash
my %hash;
while (my $line = <FH> ) {
    chomp $line;
    my ($read_name, $frags) = split(/\t/, $line);
    my @frags = split(/_/, $frags);
    my %frag_hash;                                                          #copy frags from array to frag_hash to remove duplicates
    foreach my $frag (@frags) { 
        $frag_hash{$frag} = 1;
    }
    my @frags_sorted;                                                       #sort keys in frag_hash and copy to frags_sorted array
    foreach my $X (sort keys %frag_hash) {
        push (@frags_sorted, $X);
    }
    if ($#frags_sorted == 1) {                                              #triplicates
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[1]}) {         #if new combination, add to hash with value 1
            $hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;                 #if combination exists, add 1 to count in hash               
        }
    }
    if ($#frags_sorted == 2) {                                              #quadruplicates 
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[1]}) {         
            $hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
        }
    }
    if ($#frags_sorted == 3) {                                              #quintuplicates
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[1]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[3]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[3]} += 1;
        }
        unless (exists $hash{$frags_sorted[2]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[2]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[2]}{$frags_sorted[3]} += 1;
        }
    }
    if ($#frags_sorted == 4) {                                              #sextuplicates
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[1]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[3]} += 1;
        }
        unless (exists $hash{$frags_sorted[0]}{$frags_sorted[4]}) {
            $hash{$frags_sorted[0]}{$frags_sorted[4]} = 1;
        }
        else {
            $hash{$frags_sorted[0]}{$frags_sorted[4]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[2]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[3]} += 1;
        }
        unless (exists $hash{$frags_sorted[1]}{$frags_sorted[4]}) {
            $hash{$frags_sorted[1]}{$frags_sorted[4]} = 1;
        }
        else {
            $hash{$frags_sorted[1]}{$frags_sorted[4]} += 1;
        }
        unless (exists $hash{$frags_sorted[2]}{$frags_sorted[3]}) {
            $hash{$frags_sorted[2]}{$frags_sorted[3]} = 1;
        }
        else {
            $hash{$frags_sorted[2]}{$frags_sorted[3]} += 1;
        }
        unless (exists $hash{$frags_sorted[2]}{$frags_sorted[4]}) {
            $hash{$frags_sorted[2]}{$frags_sorted[4]} = 1;
        }
        else {
            $hash{$frags_sorted[2]}{$frags_sorted[4]} += 1;
        }
        unless (exists $hash{$frags_sorted[3]}{$frags_sorted[4]}) {
            $hash{$frags_sorted[3]}{$frags_sorted[4]} = 1;
        }
        else {
            $hash{$frags_sorted[3]}{$frags_sorted[4]} += 1;
        }
    }
}

#Print to output file
foreach my $X (sort keys %hash) {
    foreach my $Y (sort keys %{$hash{$X}}) {
        print FHOUT "$X\t$Y\t$hash{$X}{$Y}\n";
    }
}


