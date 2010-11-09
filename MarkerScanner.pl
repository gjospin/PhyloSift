#!/usr/bin/perl -w


#######################################
#                                      
# Author : Guillaume Jospin                                      
# Date : November 2010
#
# Purpose : Recreates the various step from AMPHORA for updated softwares as of March 2010
#           Uses parts of the code from the AMPHORA scripts.         
#                                     
# Usage : perl Amphora_pipeline.pl                                      
#                                     
#
#######################################  






use warnings;
use strict;
use Cwd;
use Getopt::Long;


my $usage = qq~
Usage : $0 <options> <marker_name> <sequences_to_scan>

~;

my $threadNum = 1;
my $clean = 0;

GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
    ) || die $usage;

die $usage unless ($ARGV[0] && $ARGV[1]);

my $markerName = $ARGV[0];
my $sequencesFile = $ARGV[1];
my $workingDir = getcwd;


#check if the marker exists.
print "$markerName was not found\n" unless (-e "$workingDir/Markers/$markerName.faa");

#check if the temporary directory exists, if it doesn't create it.
`mkdir $workingDir/Amph_temp` unless (-e "$workingDir/Amph_temp");

#check if the sequence file exists
print "$sequencesFile was not found \n" unless (-e "$workingDir/$sequencesFile" || -e "$sequencesFile");
