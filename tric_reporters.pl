#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;

###############################################################################################################################################################################################

#This script performs the analysis of a Tri-C experiment. It takes as input the aligned sam file, a file describing the capture oligonucleotides
#targeting the viewpoints of interest and a file with the coordinates of the restriction fragments (of correspondong enzyme) in the genome.
#It removes duplicate reads and performs a proximity exclusion (all mapped restriction fragments that fall within a 1kb window from the viewpoint will
#be removed) and maps the reporters to the restriction fragments.
#See http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/UserManualforCaptureCanalysis.pdf for more information about the format of the
#required input files.

#Output:
#1. Report with basic statistics: total aligned reads, duplicates, proximity exclusion, reporters in cis and trans
#2. General txt file with all reporter fragments in the (format: read_name \t capture_oligo \t frag1-coord_frag2-coord_...)
#3. Oligo-specific txt file with reporter fragments in the (format read_name \t frag1-coord_frag2-coord_...)  
#4. Oligo-specific wig file
#5. Oligo-specific gff file

#example of run command:
#module load ucsctools
#nohup perl tric_reporters.pl -sam /t1-data1/WTSA_Dev/oudelaar/tric/pilot/analysis/Rep1_dig.sam -o /t1-data1/WTSA_Dev/oudelaar/tric/pilot/analysis/tric_oligo_file.txt -r /t1-data1/WTSA_Dev/oudelaar/scripts/mm9_nlaIII_coordinates.txt -name Rep1 -pf /public/hugheslab/ &

###############################################################################################################################################################################################

&GetOptions
(
	"sam=s" =>\ my $sam_file, 	    # -sam          sam file
    "o=s" => \my $oligo_file,       # -o            oligo file (CCanalyser format)           
    "r=s" =>\ my $dig_genome,       # -r            file with restriction coordinates in genome 
    "name=s" =>\ my $name,          # -name         name of the experiment, can be whatever you like, will be used for names of output files
    "pf=s" =>\ my $public_folder,   # -pf		    your public folder (e.g. /public/username)
);

#open filehandles
open (SAMFH, $sam_file) or die "can't open sam file $sam_file";
open (OLIFH, $oligo_file) or die "can't open oligo file $oligo_file";
open (GENFH, $dig_genome) or die "can't open genome file with restriction coordinates $dig_genome";

my $path = "undef"; 
if ($sam_file =~ /(.*\/)(\V++)/) {
    $path = $1;
    }
my $dir = "$path/$name\_tri_CC/";
unless (-d $dir) {
    mkdir $dir;
}

my $out_fragments = "$dir/$name\_all_fragments.txt";
my $report = "$dir/$name\_report.txt";

open (OUT, ">$out_fragments") or die "can't open output file $out_fragments";
open (REP, ">$report") or die "can't open output report file $report";

#store oligo coordinates in hash
my %oligos;
while (my $line = <OLIFH>) {
    chomp $line;
    my ($id, $chr1, $start1, $stop1, $chr2, $start2, $stop2, $SNP_loc, $SNP_base) = split (/\t/, $line);
    $oligos{$id}{"chr"}=$chr1;
    $oligos{$id}{"start"}=$start1;
    $oligos{$id}{"stop"}=$stop1;
    $oligos{$id}{"start_prox"}=$start2;
    $oligos{$id}{"stop_prox"}=$stop2;
}

#store restriction fragment coordinates in hash and sort in ascending order, to use for binary search
my %RE_hash;
while (my $line = <GENFH>) {
    chomp $line;
    my ($chr, $start, $stop) = split (/\W/, $line);
    push @{$RE_hash{$chr}}, $start;
}

my @chr_index = keys %RE_hash;
foreach my $chr (@chr_index) {
    @{$RE_hash{$chr}} = sort {$a <=> $b} @{$RE_hash{$chr}};
}

#store read info (name, pair-mate, fragment number [after in silico digest]) and coordinates in hash, and check if the read contains a viewpoint
my %data_hash; my %counter; 
while (my $line = <SAMFH>) {
    chomp $line;
    my ($name, $flag, $chr, $start, $map_qual, $cigar, $mate_ref, $mate_start, $mate_insert, $seq, $seq_qual, $opt1, $opt2, $opt3) = split (/\t/, $line);
    my $read_name = "undefined"; my $stop = 0;
    if ($chr =~ /chr.*/ and $name =~ /(.*):RF(\d++):(\d++)$/) {
        $read_name = $1; my $mate = $2; my $frag_nr= $3;
        $chr =~ /chr(.*)/; $chr = $1;
        $data_hash{$read_name}{$mate}{$frag_nr}{"chr"} = $chr;
        $data_hash{$read_name}{$mate}{$frag_nr}{"start"} = $start;
    if ($cigar =~/(\d++)(.*)/) {
        $stop = $1 + $start;                                                        #first number in cigar string ($1) is alignment length; if case of indels, the part of the read 
        $data_hash{$read_name}{$mate}{$frag_nr}{"stop"}= $stop;                     #until first indel is used to map the read to the restriction fragment
        my $coord = "$chr:$start-$stop";
        $data_hash{$read_name}{"coord"} .= $coord;                                  #generate coordinate string of all fragments in read for duplicate removal
        }
    foreach my $id (keys %oligos) {                                                 #check if read contains a viewpoint
        if ($chr eq $oligos{$id}{"chr"} and $start >= $oligos{$id}{"start"} - 1 and $stop <= $oligos{$id}{"stop"} + 1) {
            $data_hash{$read_name}{"label"} = $id;
            }
        }
    }
}

#remove duplicates
my %coord_hash;             
foreach my $key (keys %data_hash) {
    $counter{"1Total aligned reads"}++;
    if (exists $coord_hash {$data_hash{$key}{"coord"}}) {
        $counter{"3Duplicated aligned reads"}++;
        delete $data_hash{$key};
    }
    else {
    $coord_hash{$data_hash{$key}{"coord"}} = 1;
    $counter{"2Unique aligned reads"}++;
    }
}

#maps reads to restriction fragment using the binary_search subroutine, generate general fragment output file in format: readnr capture_oligo   frag1_frag2_... (in order),
#and generate %frag_hash for oligo specific fragment and wig output
my %frag_hash;
my $chr = 0;
foreach my $read_name (keys %data_hash) {
    foreach my $oligo (keys %oligos) {
        $chr = $oligos{$oligo}{"chr"};
        if (exists($data_hash{$read_name}{"label"}) and $data_hash{$read_name}{"label"} eq $oligo) {
            print OUT "$read_name\t$oligo\t";
            $counter{"4aUnique reads with capture total"}++;
            $counter{"4bUnique reads with capture $oligo"}++; 
            for (my $mate = 1; $mate < 3; $mate++) {                                                                                                        #sort by pair mate
                foreach my $frag (sort {$data_hash{$read_name}{$mate}{$a} <=> $data_hash{$read_name}{$mate}{$b}} keys %{$data_hash{$read_name}{$mate}}) {   #sort by frag nr
                    #viewpoint fragments
                    if ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $data_hash{$read_name}{$mate}{$frag}{"start"} >= $oligos{$oligo}{"start"} - 1 and $data_hash{$read_name}{$mate}{$frag}{"stop"} <= $oligos{$oligo}{"stop"} + 1) {
                        $data_hash{$read_name}{$mate}{$frag}{"label"}= "capture $oligo";
                    }
                    #proximity excluded fragments
                    else {
                        my ($start_frag, $end_frag) = binary_search(\@{$RE_hash{$chr}}, $data_hash{$read_name}{$mate}{$frag}{"start"}, $data_hash{$read_name}{$mate}{$frag}{"stop"}, \%counter);                    #map read to restriction fragment
                        unless ($start_frag =~ /error/) {
                            if ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag <= $oligos{$oligo}{"start_prox"} - 1 and $end_frag >= $oligos{$oligo}{"start_prox"} - 1) {    #proximity exclusion if end of fragment is in proximity zone 
                                $data_hash{$read_name}{$mate}{$frag}{"label"} = "prox_excl $oligo";
                                $counter{"5aProximity excluded fragments total"}++;
                                $counter{"5bProximity excluded fragments $oligo"}++;
                            }
                            elsif ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag >= $oligos{$oligo}{"start_prox"} - 1 and $end_frag <= $oligos{$oligo}{"stop_prox"} + 1) {  #proximity exclusion if entire fragment is in proximity zone
                                $data_hash{$read_name}{$mate}{$frag}{"label"} = "prox_excl $oligo";
                                $counter{"5aProximity excluded fragments total"}++;
                                $counter{"5bProximity excluded fragments $oligo"}++;
                            }
                            elsif ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag <= $oligos{$oligo}{"stop_prox"} - 1 and $end_frag >= $oligos{$oligo}{"stop_prox"} + 1) {   #proximity exclusion if start of fragment is in proximity zone
                                $data_hash{$read_name}{$mate}{$frag}{"label"} = "prox_excl $oligo";
                                $counter{"5aProximity excluded fragments total"}++;
                                $counter{"5bProximity excluded fragments $oligo"}++;
                            }
                            #reporter fragments (remaining fragments)
                            else {
                                $data_hash{$read_name}{$mate}{$frag}{"label"}= "reporter $oligo";
                                print OUT "$data_hash{$read_name}{$mate}{$frag}{\"chr\"}:$start_frag-$end_frag\_";
                                push (@{$frag_hash{$oligo}{$read_name}}, "$data_hash{$read_name}{$mate}{$frag}{\"chr\"}:$start_frag-$end_frag");
                                $counter{"6aReporter fragments total"}++;
                                $counter{"6bReporter fragments $oligo"}++;
                                if ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"}) {
                                    $counter{"6cReporter fragments $oligo cis"}++;
                                }
                                else {
                                    $counter{"6dReporter fragments $oligo trans"}++;
                                }
                            }
                        }
                    }
                }
            }
            print OUT "\n";
        }
    }
}

#print basic statistics to report file
foreach my $key (sort keys %counter) {
    printf REP "%-8s %s\n", $key, $counter{$key};
}

print REP "\nTracks:\n";

#generate oligo-specific fragment output file 
foreach my $oligo (keys %frag_hash) {
    my $frag_output = "$dir/$name\_$oligo\_fragments.txt";
    open (FRAG_OUT, ">$frag_output") or die "can't open fragment output file $frag_output";
    foreach my $read (keys %{$frag_hash{$oligo}}) {
        print FRAG_OUT "$read\t";
        my $length = scalar @{$frag_hash{$oligo}{$read}};
        for (my $i = 0; $i < $length; $i++ ) {
            my ($chr, $start, $stop) = split(/\W/, ${frag_hash{$oligo}{$read}}[$i]);            
            print FRAG_OUT "$chr:$start-$stop\_";
        }
        print FRAG_OUT "\n";   
    }
}

#generate gff hash 
my %gff_hash;
foreach my $oligo (keys %frag_hash) {
    foreach my $read (keys %{$frag_hash{$oligo}}) {
        my $length = @{$frag_hash{$oligo}{$read}};
        for(my $i = 0; $i < $length; $i++) {
            my $frag = ${frag_hash{$oligo}{$read}}[$i];
            ${$gff_hash{$oligo}{$frag}}++;            
        }
    }
}

#generate oligo-specific gff output file 
foreach my $oligo (keys %gff_hash) {
    my $gff_output = "$dir/$name\_$oligo\.gff";
    open (GFF_OUT, ">$gff_output") or die "can't open fragment output file $gff_output";
    foreach my $coord (sort keys %{$gff_hash{$oligo}}) {
        my ($chr, $dpn_start, $dpn_stop) = split(/\W/, $coord);
        $dpn_start = $dpn_start + 2;
        $dpn_stop = $dpn_stop + 2;
        unless ($chr =~ /M|Y/) {                                                                            #exclude reporters on chrM and chrY: binary search gives unreported error
            print GFF_OUT "chr$chr\tco-cap\t$name\t$dpn_start\t$dpn_stop\t${$gff_hash{$oligo}{$coord}}\t+\t0\t.\n";
        }
    }
}

#generate tracks
frag_to_wig(\%frag_hash,"");

###############################################################################################################################################################################################

sub binary_search {
    my ($chr_array, $start, $stop, $counter_hash) = @_;
    my $mid_value = ($start + $stop)/2;
    my $first_pos = 0;
    my $last_pos = scalar @$chr_array - 1; 
    my $counter =0;
    if (($mid_value < $$chr_array[$first_pos]) or ($mid_value > $$chr_array[$last_pos])) {
        $$counter_hash{"Binary search error: search outside range of restriction enzyme coords"}++;
        return ("error", "error")}
    for (my $i=0; $i<99; $i++) {
        my $mid_search = int(($first_pos + $last_pos) / 2);
        if ($$chr_array[$mid_search] > $$chr_array[$mid_search+1]) {
            $$counter_hash{"Binary search error: restriction enzyme array coordinates not in ascending order"}++;
            return ("error", "error")}
        if (($$chr_array[$mid_search] <= $mid_value) and ($$chr_array[$mid_search+1] > $mid_value)) {    #maps the mid point of the read to a fragment
            if (($$chr_array[$mid_search] <= $start+2) and ($$chr_array[$mid_search+1] >= $stop-2)) {    #checks the whole read is on the fragment +/-2 to allow for the dpnII overlaps
                return ($$chr_array[$mid_search], $$chr_array[$mid_search+1]-1)
                }
            else {
                $$counter_hash{"Binary search error: fragment overlaps multiple restriction sites"}++;
                return ("error", "error");
                }
        }       
        elsif ($$chr_array[$mid_search] > $mid_value) {
            $last_pos = $mid_search-1;
            }    
        elsif ($$chr_array[$mid_search] < $mid_value) {
            $first_pos = $mid_search+1;
            }
        else {
            $$counter_hash{"Binary search error: end of loop reached"}++;
            }
    }
    $$counter_hash{"Binary search error: couldn't map read to fragments"}++;
    return ("error", "error")
}

sub frag_to_wig {                                                                       #hashformat {id}{read}{1:1000-2000}
my ($hashref, $file_name) = @_; 
    foreach my $oligo (keys %$hashref) {
    foreach my $read (keys %{$$hashref{$oligo}}) {
        my $length = scalar @{$$hashref{$oligo}{$read}};
        for (my $i = 0; $i < $length; $i++ ) {
            my ($chr, $start, $stop) = split(/\W/, $$hashref{$oligo}{$read}[$i]);
            #$chr =~ s/chr//gi;
            if ($chr eq $oligos{$oligo}{"chr"}) {                                       #only plot in cis to make script run faster
                my @range = ($start..$stop);
                for(my $j = 0; $j < $#range; $j++) {
                    ${$coord_hash{$oligo}{$chr}{$range[$j]}}++;
                    }
            }
            }
        }
    my $tracks_out = "$name\_$oligo$file_name";
    open (TRACKS_OUT, ">$dir/$tracks_out.wig");
    foreach my $chr (sort {$coord_hash{$oligo}{$a} <=> $coord_hash{$oligo}{$b}} keys %{$coord_hash{$oligo}}) {
        print TRACKS_OUT "variableStep  chrom=chr$chr\n";
        foreach my $coord (sort {$a <=> $b} keys %{$coord_hash{$oligo}{$chr}}) {
            my $count = ${$coord_hash{$oligo}{$chr}{$coord}};
            print TRACKS_OUT "$coord\t$count\n";
        }
    }
    system ("wigToBigWig -clip $dir/$tracks_out.wig /t1-data/user/config/bigwig/mm9_sizes.txt $dir/$tracks_out.bw") == 0 or die "couldn't bigwig files\n";
    system ("mv $dir/$tracks_out.bw $public_folder") == 0 or die "couldn't move files\n";		
    system ("chmod 755 $public_folder/$tracks_out.bw") == 0 or die "couldn't chmod files\n";   
    print REP "track type=bigWig name=\"$tracks_out\" description=\"c-trap $tracks_out\" bigDataUrl=http://sara.molbiol.ox.ac.uk$public_folder/$tracks_out.bw\n";
    }
}
    

