#!/usr/bin/perl
package Phylosift::Phylosift;
use 5.006;
use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use File::Basename;
use Carp;
use Phylosift::Utilities qw(:all);
use Phylosift::MarkerAlign;
use Phylosift::pplacer;
use Phylosift::Summarize;
use Phylosift::FastSearch;
use Phylosift::Benchmark;
use Phylosift::BeastInterface;
use Phylosift::Comparison;
use Phylosift::MarkerBuild;

=head2 new

    Returns : Phylosift project object
    Args : pair,readsFile(,readsFile_2);

=cut

sub new {
	my $self = {};
	$self->{"fileName"}    = undef;
	$self->{"workingDir"}  = undef;
	$self->{"mode"}        = undef;
	$self->{"readsFile"}   = undef;
	$self->{"readsFile_2"} = undef;
	$self->{"tempDir"}     = undef;
	$self->{"fileDir"}     = undef;
	$self->{"blastDir"}    = undef;
	$self->{"alignDir"}    = undef;
	$self->{"treeDir"}     = undef;
	$self->{"dna"}         = undef;
	$self->{"updated"}     = 0;
	$self->{"coverage"}    = undef;
	$self->{"isolate"}     = 0;
	$self->{"threads"}     = 1;
	bless($self);
	return $self;
}

=head2 initialize
    
    Initializes the variables for the Phylosift object
    Using the standard pathnames and the filename

=cut

sub initialize {
	my $self      = shift;
	my $mode      = shift;
	my $readsFile = shift;
	debug "READSFILE\t" . $readsFile . "\n";
	my $readsFile_2 = "";
	if ( scalar(@_) == 1 ) {
		debug "FOUND a second file\n";
		$readsFile_2 = shift;
	}
	my $position = rindex( $readsFile, "/" );
	$self->{"fileName"}    = substr( $readsFile, $position + 1, length($readsFile) - $position - 1 );
	$self->{"workingDir"}  = getcwd;
	$self->{"mode"}        = $mode;
	$self->{"readsFile"}   = $readsFile;
	$self->{"readsFile_2"} = $readsFile_2;
	$self->{"tempDir"}     = $self->{"workingDir"} . "/PS_temp";
	$self->{"fileDir"}     = $self->{"tempDir"} . "/" . $self->{"fileName"};
	$self->{"blastDir"}    = $self->{"fileDir"} . "/blastDir";
	$self->{"alignDir"}    = $self->{"fileDir"} . "/alignDir";
	$self->{"treeDir"}     = $self->{"fileDir"} . "/treeDir";
	$self->{"dna"}         = 0;
	$self->{"updated"}     = 0;
	$self->{"16s"}         = 0;
	return $self;
}

=head2 getReadsFile

    returns the file name for the reads

=cut

sub getReadsFile {
	my $self = shift;
	return $self->{"readsFile"};
}

=head1 NAME

Phylosift::Phylosift - Implements core functionality for Phylosift

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Phylosift::Phylosift;

    my $foo = Phylosift::Phylosift->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 run
    
    Args : $force - if used at hte same time as the all mode it removes the temp directory if it already exists
           $custom - reads in marker names from a custom list. if "" then use all markers in the marker directory
           $continue - if the mode != all and continue = 1 then finish the pipeline otherwise only execute the step specified
           $isolateMode - allows sequences to hit multiple markers (used when running isolate genomes)

    Runs the PhyloSift pipeline according to the functions passed as arguments

=cut

my $continue = 0;
my ( $mode, $readsFile, $readsFile_2, $fileName, $tempDir, $fileDir, $blastDir, $alignDir, $treeDir ) = "";
my ( $sec,  $min,       $hour,        $mday,     $mon,     $year,    $wday,     $yday,     $isdst )   = 0;
my $workingDir = getcwd;

#where everything will be written when PhyloSift is running
sub run {
	my $self     = shift;
	my $force    = shift;
	my $custom   = shift;
	my $continue = shift;
	debug "force : $force\n";
	Phylosift::Utilities::print_citations();
	start_timer("START");
	$self->readPhylosiftConfig();
	$self->runProgCheck();
	Phylosift::Utilities::data_checks(self=>$self);
	$self->fileCheck();
	$self->directoryPrep($force);
	$self->{"readsFile"} = $self->prepIsolateFiles( $self->{"readsFile"} ) if $self->{"isolate"} == 1;

	#create a file with a list of markers called markers.list
	debug "CUSTOM = " . $custom . "\n";
	my @markers = Phylosift::Utilities::gather_markers( self=>$self, marker_file => $custom );
	if($self->{"extended"}){
		@markers = Phylosift::Utilities::gather_markers( self=>$self, path => $Phylosift::Utilities::markers_extended_dir );
	}
	debug "@markers\n";
	debug "MODE :: " . $self->{"mode"} . "\n";
	if ( $self->{"mode"} eq 'search' || $self->{"mode"} eq 'all' ) {
		$self = $self->runSearch( $continue, $custom, \@markers );
		debug "MODE :: " . $self->{"mode"} . "\n";
	}
	debug "MODE :: " . $self->{"mode"} . "\n";
	if ( $self->{"mode"} eq 'align' || $self->{"mode"} eq 'all' ) {
		$self = $self->runMarkerAlign( $continue, \@markers );
	}
	if ( $self->{"mode"} eq 'placer' || $self->{"mode"} eq 'all' ) {
		croak "No marker gene hits found in the input data. Unable to reconstruct phylogeny and taxonomy." if ( scalar(@markers) == 0 );
		$self = $self->runPplacer( $continue, \@markers );
	}
	if ( $self->{"mode"} eq 'summary' || $self->{"mode"} eq 'all' ) {
		$self = $self->taxonomyAssignments( $continue, \@markers );
	}
	if ( $self->{"mode"} eq 'benchmark' ) {
		$self = $self->benchmark();
	}
	if ( $self->{"mode"} eq 'compare' ) {
		$self = $self->compare();
	}
	if ( $self->{"mode"} eq 'index' ) {
		Phylosift::Utilities::index_marker_db( self=>$self, markers=>\@markers, path=>$Phylosift::Utilities::marker_dir );
		my @extended_markers = Phylosift::Utilities::gather_markers( self=>$self, path => $Phylosift::Utilities::markers_extended_dir );
		Phylosift::Utilities::index_marker_db( self=>$self, markers=>\@extended_markers, path=>$Phylosift::Utilities::markers_extended_dir );
	}
	if( $self->{"mode"} eq 'build_marker'){
		Phylosift::MarkerBuild::build_marker(alignment=>$ARGV[1], name=>$ARGV[2], cutoff=>$ARGV[3]);
	}
}

=head2 function2

    Reads the Phylosift configuration file and assigns the file paths to the required directories

=cut

# reads the Phylosift configuration file
sub readPhylosiftConfig {
	my $self          = shift;
	my $custom_config = $self->{"configuration"};

	# first get the install prefix of this script
	my $scriptpath = dirname($0);

	# try first a config in the script dir, in case we're running from
	# a dev directory.  then config in system dir, then user's home.
	# let each one override its predecessor.
	{

		package Phylosift::Settings;
		do "$scriptpath/phylosiftrc";
		do "$scriptpath/../phylosiftrc";
		do "$scriptpath/../etc/phylosiftrc";
		do "$ENV{HOME}/.phylosiftrc";
		do $custom_config if defined $custom_config;
	}
	return $self;
}

=head2 fileCheck

    Checks if the files passed to the Phylosift object exist and are not empty

=cut

sub fileCheck {
	my $self = shift;
	if ( !-e $self->{"readsFile"} ) {
		die $self->{"readsFile"} . "  was not found \n";
	}

	#check if the input file is a file and not a directory
	if ( !-f $self->{"readsFile"} ) {
		die $self->{"readsFile"} . " is not a plain file, could be a directory\n";
	}
	if ( $self->{"readsFile_2"} ne "" ) {

		#check the input file exists
		if ( !-e $self->{"readsFile_2"} ) {
			die $self->{"readsFile_2"} . " was not found\n";
		}

		#check if the input file is a file and not a directory
		if ( !-e $self->{"readsFile_2"} ) {
			die $self->{"readsFile_2"} . " is not a plain file, could be a directory\n";
		}
	}
	return $self;
}

=head2 runProgCheck

    Runs a check on the programs that will be used through the pipeline to make sure they are
    available to the user and are the versions Phylosift was tested with.

=cut

sub runProgCheck {
	my $self = shift;

	#check if the various programs used in this pipeline are installed on the machine
	my $progCheck = Phylosift::Utilities::programChecks($self);
	if ( $progCheck != 0 ) {
		croak "A required program was not found during the checks aborting\n";
	} elsif ( $progCheck == 0 ) {
		debug "All systems are good to go, continuing the screening\n";
	}
	return $self;
}

=head2 prepIsolateFiles

=item *

Process an input file for isolate mode.  Creates a temporary input file that contains only a single
sequence entry.  TODO: this could be made more elegant by storing a mapping between sequence entries
and isolate names in memory and using that at later stages of the pipeline

=back

=cut

sub prepIsolateFiles {
	my $self = shift;
	open( OUTFILE, ">" . $self->{"fileDir"} . "/isolates.fasta" );
	while ( my $file = shift ) {
		open( ISOLATEFILE, $file ) || croak("Unable to read $file\n");
		my $fname = $self->{"fileDir"} . "/" . basename($file);
		debug "Operating on isolate file $fname\n";
		print OUTFILE ">" . basename($file) . "\n";
		while ( my $line = <ISOLATEFILE> ) {
			next if $line =~ /^>/;
			print OUTFILE $line;
		}
		close ISOLATEFILE;
	}
	close OUTFILE;
	$self->{"readsFile"} = "isolates.fasta";
	return $self->{"fileDir"} . "/isolates.fasta";
}


=head2 directoryPrep

    Prepares the temporary Phylosift directory by deleting old runs and/or creating the correct directory structure
    
=cut

sub directoryPrep {
	my $self  = shift;
	my $force = shift;

	#    print "FORCE DIRPREP   $force\t mode   ".$self->{"mode"}."\n";
	#    exit;
	#remove the directory from a previous run
	if ( $force && $self->{"mode"} eq 'all' ) {
		debug( "deleting an old run\n", 0 );
		my $dir = $self->{"fileDir"};
		`rm -rf $dir`;
	} elsif ( -e $self->{"fileDir"} && $self->{"mode"} eq 'all' ) {
		croak(   "A previous run was found using the same file name aborting the current run\n"
			   . "Either delete that run from "
			   . $self->{"fileDir"}
			   . ", or force overwrite with the -f command-line option\n" );
	}

	#check if the temporary directory exists, if it doesn't create it.
	`mkdir $self->{"tempDir"}` unless ( -e $self->{"tempDir"} );

	#create a directory for the Reads file being processed.
	`mkdir $self->{"fileDir"}`  unless ( -e $self->{"fileDir"} );
	`mkdir $self->{"blastDir"}` unless ( -e $self->{"blastDir"} );
	`mkdir $self->{"alignDir"}` unless ( -e $self->{"alignDir"} );
	`mkdir $self->{"treeDir"}`  unless ( -e $self->{"treeDir"} );
	return $self;
}

=head2 taxonomyAssignments

    performs the appropriate checks before running the taxonomy classification parts of the pipeline

=cut

sub taxonomyAssignments {
	my $self        = shift;
	my $continue    = shift;
	my $markListRef = shift;
	Phylosift::Utilities::start_timer("taxonomy assignments");
	Phylosift::Summarize::summarize( $self, $markListRef );
	Phylosift::Utilities::end_timer("taxonomy assignments");
	return $self;
}

=head2 benchmark

=cut

sub benchmark {
	my $self = shift;
	debug "RUNNING Benchmark\n";
	Phylosift::Benchmark::runBenchmark( $self, "./" );
}

=head2 compare

=cut

sub compare {
	my $self = shift;
	debug "RUNNING Compare\n";
	Phylosift::Comparison::compare( $self, "./" );
}

=head2 runPplacer

    Performs the appropriate checks before runing the Read placement parts of the pipeline

    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub runPplacer {
	my $self        = shift;
	my $continue    = shift;
	my $markListRef = shift;
	debug "PPLACER MARKS @{$markListRef}\n";
	Phylosift::Utilities::start_timer("runPPlacer");
	Phylosift::pplacer::pplacer( $self, $markListRef );
	Phylosift::Utilities::end_timer("runPPlacer");

	if ( $continue != 0 ) {
		$self->{"mode"} = 'summary';
	}
	return $self;
}

=head2 runMarkerAlign

    Run the Hmm hit verification and the hit alignments to prep the hits for Pplacer
    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub runMarkerAlign {
	my $self     = shift;
	my $continue = shift;
	my $markRef  = shift;
	Phylosift::Utilities::start_timer("Alignments");

	#clearing the alignment directory if needed
	my $alignDir = $self->{"alignDir"};
	`rm $alignDir/*` if (<$alignDir/*>);

	#Align Markers
	my $threadNum = 1;
	Phylosift::MarkerAlign::MarkerAlign( $self, $markRef );

	#    Phylosift::BeastInterface::Export($self, $markRef, $self->{"fileDir"}."/beast.xml");
	Phylosift::Utilities::end_timer("Alignments");
	if ( $continue != 0 ) {
		$self->{"mode"} = 'placer';
	}
	return $self;
}

=head2 runBlast

    runBlast( markerArray, readsFile)
    
    Runs reads marker classifiation of the pipeline using Blast

    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub runSearch {
	my $self          = shift;
	my $continue      = shift;
	my $custom        = shift;
	my $type          = shift;
	my $markerListRef = shift;
	Phylosift::Utilities::start_timer("runBlast");

	#clearing the blast directory
	my $blastDir = $self->{"blastDir"};
	`rm $self->{"blastDir"}/*` if (<$blastDir/*>);

	#run Blast
	Phylosift::FastSearch::RunSearch( $self, $custom, $type, $markerListRef );
	Phylosift::Utilities::end_timer("runBlast");
	if ( $continue != 0 ) {
		$self->{"mode"} = 'align';
	}
	return $self;
}

=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::Phylosift


You can also look for information at:


=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Phylosift-Phylosift>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Phylosift-Phylosift>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Phylosift-Phylosift>

=item * Search CPAN

L<http://search.cpan.org/dist/Phylosift-Phylosift/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation.

See http://dev.perl.org/licenses/ for more information.


=cut

1;    # End of Phylosift::Phylosift
