#!/usr/bin/perl -w
use strict;
use warnings;
######################################
#
# Author : Aaron Darling and Guillaume Jospin
#
#
#
#
# demultiplex_dualBC.pl <illumina_directory> <barcode_list> <samples_mapping> <output_dir> <output_file_core_name>
#
# Illumina_directory : Directory where the illumina files are located for the dual barcode demultiplex
#                      For a MiSeq run, files should look like :  
#                      XTDM1_NoIndex_L001_R1_001.fastq.gz
#                      XTDM1_NoIndex_L001_R2_001.fastq.gz
#                      XTDM1_NoIndex_L001_R3_001.fastq.gz
#                      XTDM1_NoIndex_L001_R4_001.fastq.gz
# Barcode_list : list of the barcode name and barcode.  <barcode_label><tab><barcode>
# samples_mapping : List of the samples names to barcode combination.  <barcode1>.<barcode2><tab><name_to_use>
#                   Any spaces in the name to use column will be replaced by a _
#
# output_file_core_name : Core name for all the samples from the sequencing run.
# output_dir : Directory to output files to. Need a CLEAN directory. 
#              Files will be appended so files from an old run need to be removed
#
# Printing to STDERR the read counts for each barcode
# Mismatched_barcode : one of the two mates had a barcode that was not recognized from the list
# BC1.BC2 : Whenever a barcode pair was not listed in the name mapping file
#
######################################


my $usage = "Wrong number of arguments\nUsage:\ndemultiplex_dualBC.pl <illumina_directory> <barcode_list> <samples_mapping> <output_dir> <output_file_core_name>\n";
die ("$usage") if @ARGV != 5;


#reading barcodes
my %barcode = ();
my %barcode_rev = ();

print STDERR "Reading barcodes\n";
open(INBC,$ARGV[1]);
while(<INBC>){
    chomp($_);
    my @line = split(/\t/,$_);
    $barcode{$line[1]}=$line[0];
}
close(INBC);

#reading barcodes to names mapping
my %mapping = ();
print STDERR "Reading name mapping\n";
open(INMAP,$ARGV[2]);
while(<INMAP>){
    chomp($_);
    my @line = split(/\t/,$_);
    $line[1] =~ s/\s/_/g;
    $mapping{$line[0]}=$line[1];
}
close(INBC);

my %bc_count = ();
my @files = <$ARGV[0]/*_R1_*.fastq.gz>;

foreach my $file(@files){
    print STDERR "Processing $file\n";
    $file =~ m/^(\S+)_\S\S_(\d+).fastq.gz/;
    my $core = $1;
    my $index = $2;
    open(READ1,"zcat $file |");
    open(READ2,"zcat $core"."_R4_$index.fastq.gz |");
    open(INDEX1,"zcat $core"."_R2_$index.fastq.gz |");
    open(INDEX2,"zcat $core"."_R3_$index.fastq.gz |");
    my @read1 = ();
    my @read2 = ();
    my @index1 = ();
    my @index2 = ();
    while(1){
	for(my $i=0; $i<4; $i++){
	    $read1[$i] = <READ1>;
	    $read2[$i] = <READ2>;
	    $index1[$i] = <INDEX1>;
	    $index2[$i] = <INDEX2>;
	}
	last if !defined($read1[0]);
	my $BC1 = "";
	my $BC2 = "";
	my $i1 = $index1[1];
	my $i2 = $index2[1];
	chomp($i1);
	chomp($i2);
	if(exists $barcode{$i1}){
	    $BC1 = $barcode{$i1};
	}
	if(exists $barcode{$i2}){
	    $BC2 = $barcode{$i2};
        }
	if($BC1 eq "" || $BC2 eq ""){
	    $bc_count{mismatch} = 1 if ! exists $bc_count{mismatch};
	    $bc_count{mismatch} ++ if exists $bc_count{mismatch};
	    open(OUT1,">>$ARGV[3]/$ARGV[4]"."_mismatched_barcodes"."_1.fastq");
	    open(OUT2,">>$ARGV[3]/$ARGV[4]"."_mismatched_barcodes"."_2.fastq");
	    print OUT1 @read1;
	    print OUT2 @read2;
	    close(OUT1);
	    close(OUT2);
	}else{
	    my $code = $BC1.".".$BC2;
	    $code = $mapping{"$BC1.$BC2"} if exists ($mapping{"$BC1.$BC2"});
	    $bc_count{$code}=1 if ! exists $bc_count{$code};
            $bc_count{$code} ++ if exists $bc_count{$code};
	    open(OUT1,">>$ARGV[3]/$ARGV[4]"."_$code"."_1.fastq");
	    open(OUT2,">>$ARGV[3]/$ARGV[4]"."_$code"."_2.fastq");
	    print OUT1 @read1;
	    print OUT2 @read2;
	    close(OUT1);
	    close(OUT2);
	}
    }
}

foreach my $code (sort {$bc_count{$a} cmp $bc_count{$b} } keys %bc_count){
    print STDERR $code."\t".$bc_count{$code}."\n";
}
