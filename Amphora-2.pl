#!/usr/bin/perl -w

use warnings;
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;


my $usage = qq~
Usage: $0 <options> <reads_file>

~;

my $threadNum = 1;
my $clean = 0;
my $euk = 0;
my $arc = 0;
my $bac = 0;
my $custom = "";



GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
	   "euk" => \$euk,
	   "bac" => \$bac,
	   "arc" => \$arc,
	   "custom=s" => \$custom, #need a file containing the marker names to use without extensions ** marker names shouldn't contain '_'
    ) || die $usage;

#euk, bac and arc represent which markers to run depending on their domain - run all by default
#custom overrides the 3 previous options (display warning if both are present OR die and require correct input)

#check for an input file on the command line
die $usage unless ($ARGV[0]);

my $workingDir = getcwd;
my $readsFile = $ARGV[0];
my %seq =();
#check the input file exists
print "$readsFile was not found \n" unless (-e "$workingDir/$readsFile" || -e "$readsFile");
#my $fileName = "";
my $position = rindex($readsFile,"/");
my $fileName = substr($readsFile,$position+1,length($readsFile)-$position-1);
my $tempDir = "$workingDir/Amph_temp";
my $fileDir = "$tempDir/$fileName";
#check if the temporary directory exists, if it doesn't create it.
`mkdir $tempDir` unless (-e "$tempDir");

#create a directory for the Reads file being processed.
`mkdir $fileDir` unless (-e "$fileDir");
#clear the directory of the fileName already exists
`rm -r $fileDir/*` unless (-e "$fileDir/*");


#make a blastable DB for all markers (If the blast time is too long, pick representatives, possible to add as an option)


#create a file with a list of markers
my @markers = ();

if($custom ne ""){
    #gather a custom list of makers
    `cp $custom $fileDir/markers.list`;
    
#    `perl $workingDir/run_blast.pl $custom $readsFile`;


}else{
    #gather all markers
    #LATER : add differentiation for euk - bac - arc
    open(markersOUT,">$fileDir/markers.list");
    my @files = <$workingDir/markers/*.faa>;
    foreach my $file (@files){
	$file =~ m/\/(\w+).faa/;
	push(@markers,$1);
	print markersOUT "$1\n";
    }
    close(markersOUT);
    
}


#run Blast

`perl $workingDir/run_blast.pl --threaded=$threadNum $fileDir/markers.list $readsFile`;

#Align Markers
`perl $workingDir/MarkerAlign.pl --threaded=$threadNum $fileDir/markers.list $readsFile`;

# Run Pplacer
`perl $workingDir/Run_Pplacer.pl --threaded=$threadNum $fileDir/markers.list $readsFile`


#TODO : 

# set  up alternate blast parameters if dealing with short reads (add an option or by default check the size)
# add in a 6 frame translation module (ask aaron) if nucleotides are supplied
# set up check points in case the pipeline stops
# allow for a parallel use
