package Amphora2::Amphora2;

use 5.006;
use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use File::Basename;
use Amphora2::Utilities;
use Amphora2::MarkerAlign;
use Amphora2::blast qw(RunBlast);
use Amphora2::pplacer;
use Amphora2::Summarize;

=head1 NAME

Amphora2::Amphora2 - Implements core functionality for Amphora2

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Amphora2::Amphora2;

    my $foo = Amphora2::Amphora2->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 run

=cut

sub run {

	#message to output if the script isn't called properly
	my $usage = qq~
	Usage: $0 <options> <reads_file>

	~;
	my $usage2 = qq~
	Usage: $0 <options> -paired <reads_file_1> <reads_file_2>
	~;
	#euk,arc,bac,clean don't do anything at this time
	my $threadNum = 1;
	my $clean = 0;
	my $euk = 0;
	my $arc = 0;
	my $bac = 0;
	my $custom = "";
	my $force=0;
	my $pair =0;

	GetOptions("threaded=i" => \$threadNum,
		   "clean" => \$clean,
		   "euk" => \$euk,
		   "bac" => \$bac,
		   "arc" => \$arc,
		   "paired" => \$pair, # used for paired fastQ input split in 2 different files
		   "custom=s" => \$custom, #need a file containing the marker names to use without extensions ** marker names shouldn't contain '_'
		   "f" => \$force, #overrides a previous run otherwise stop
	    ) || die $usage;

	readAmphora2Config();

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR  "START : %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

	#euk, bac and arc represent which markers to run depending on their domain - run all by default
	#custom overrides the 3 previous options (currently  only read the custom marker list)

	#check for an input file on the command line
	my $readsFile_2="";
	if($pair ==0){
	    die $usage unless ($ARGV[0]);
	}else{
	    die $usage2 unless ($ARGV[0] && $ARGV[1]);
	    $readsFile_2= $ARGV[1];
	}
	my $workingDir = getcwd;
	my $readsFile = $ARGV[0];




	#check if the various programs used in this pipeline are installed on the machine
	my $progCheck = Amphora2::Utilities::programChecks();
	if($progCheck!=0){
	    print STDERR "A required program was not found during the checks aborting\n";
	    exit();
	}elsif($progCheck==0){
	    print STDERR "All systems are good to go, continuing the screening\n";
	}


	#check the input file exists
	if(!-e "$workingDir/$readsFile" && !-e "$readsFile"){
	    die "$readsFile was not found \n";
	}
	#check if the input file is a file and not a directory
	if(!-f "$workingDir/$readsFile" && !-f "$readsFile"){
	    die "$readsFile is not a plain file, could be a directory\n";
	}
	if($pair !=0){
	#check the input file exists
	    if(!-e "$workingDir/$readsFile_2" || !-e "$readsFile_2"){
		die "$readsFile_2 was not found\n";
	    }
	#check if the input file is a file and not a directory
	    if(!-f "$workingDir/$readsFile_2" || !-e "$readsFile_2"){
		die "$readsFile_2 is not a plain file, could be a directory\n";
	    }
	}


	#die "$readsFile was not found \n" unless (-e "$workingDir/$readsFile" || -e "$readsFile");
	#get the filename in case a filepath is included
	my $position = rindex($readsFile,"/");
	my $fileName = substr($readsFile,$position+1,length($readsFile)-$position-1);

	#where everything will be written when Amphora-2 is running
	my $tempDir = "$workingDir/Amph_temp";
	#where everything will be written when running Amphora-2 using this file
	my $fileDir = "$tempDir/$fileName";

	#remove the directory from a previous run
	if($force){
	    print STDERR "deleting an old run\n";
	    `rm -rf $fileDir`;
	}elsif(-e "$fileDir"){
	    print STDERR "A previous run was found using the same file name aborting the current run\n";
	    print STDERR "Either delete that run from $fileDir, or force overwrite with the -f command-line option\n";
	    exit;
	}

	#check if the temporary directory exists, if it doesn't create it.
	`mkdir $tempDir` unless (-e "$tempDir");
	#create a directory for the Reads file being processed.
	`mkdir $fileDir` unless (-e "$fileDir");
	#clear the directory of the fileName already exists
	`rm -r $fileDir/*` unless (!-e "$fileDir/*");


	#create a file with a list of markers called markers.list
	my @markers = ();

	if($custom ne ""){
	    #gather a custom list of makers
	    `cp $custom $fileDir/markers.list`;
	    
	}else{
	    #gather all markers
	    #LATER (maybe) : add differentiation for euk - bac - arc
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
	if($pair == 0){
		Amphora2::blast::RunBlast( ("--threaded=$threadNum", "$fileDir/markers.list", "$readsFile" ) );
#	    `run_blast.pl --threaded=$threadNum $fileDir/markers.list $readsFile`;
	}elsif($pair == 1){
		Amphora2::blast::RunBlast();
	    print STDERR "Starting run_blast.pl with 2 fastQ files\n";
	    `run_blast.pl --threaded=$threadNum -paired $fileDir/markers.list $readsFile $readsFile_2`;
	}

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR "Before Alignments for Markers %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
	#Align Markers
	Amphora2::MarkerAlign::MarkerAlign( ("--threaded=$threadNum", "$fileDir/markers.list", "$readsFile" ) );
#	`MarkerAlign.pl --threaded=$threadNum $fileDir/markers.list $readsFile`;
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR "After Alignments %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;


	# Run Pplacer
	Amphora2::MarkerAlign::pplacer( ("--threaded=$threadNum", "$fileDir/markers.list", "$readsFile") );
#	`Run_Pplacer.pl --threaded=$threadNum $fileDir/markers.list $readsFile`;
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR "After PPlacer %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

	# Taxonomy assignemnts
	`/home/gjospin/Amphora-2/printneighbor.R $fileDir/trees/*.num.tre > $fileDir/neighbortaxa.txt`;
	Amphora2::Summarize::summarize();
#	`/home/gjospin/Amphora-2/summarize.pl $fileDir/neighbortaxa.txt > $fileDir/taxasummary.txt`;
	#TODO : 


	# transform the .place files into tree files.
	# set  up alternate blast parameters if dealing with short reads (add an option or by default check the size)
	# set up check points in case the pipeline stops
	# add a summary file (number of hits per markers, total hits, total reads, Nucl Vs AA, time stamps for the various steps)
	# allow for a parallel use (done for Pplacer using ForkManager, other parts are fast enough that it may not be needed)
	# look into write my own blast parser and use -m6 instead of -m0 to save on temporary file sizes.  Had issues with bioperl and using the -m6 or -m7 formats from the newest blast.

}

=head2 function2

=cut

# reads the Amphora2 configuration file
sub readAmphora2Config {
	# first get the install prefix of this script
	my $scriptpath = dirname($0);
	# try first a config in the script dir, in case we're running from
	# a dev directory.  then config in system dir, then user's home.
	# let each one override its predecessor.
	{ package Amphora2Settings; do "$scriptpath/amphora2rc"; do "$scriptpath/../etc/amphora2rc"; do "$ENV{HOME}/.amphora2rc" }
}


=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-amphora2-amphora2 at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Amphora2-Amphora2>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Amphora2::Amphora2


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Amphora2-Amphora2>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Amphora2-Amphora2>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Amphora2-Amphora2>

=item * Search CPAN

L<http://search.cpan.org/dist/Amphora2-Amphora2/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of Amphora2::Amphora2
