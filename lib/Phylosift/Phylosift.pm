#!/usr/bin/env perl
package Phylosift::Phylosift;
use 5.006;
use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::NFSLock qw(uncache);
use Fcntl qw(LOCK_EX LOCK_SH);
use Carp;
use Phylosift::Utilities qw(:all);
use Phylosift::MarkerAlign;
use Phylosift::pplacer;
use Phylosift::Summarize;
use Phylosift::FastSearch;
use Phylosift::Settings;
use Phylosift::Benchmark;
use Phylosift::BeastInterface;
use Phylosift::Comparison;
use Phylosift::MarkerBuild;
use Phylosift::Simulations;

our $VERSION = "v1.0.1";

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
	$self->{"blastDir"}    = undef;
	$self->{"alignDir"}    = undef;
	$self->{"treeDir"}     = undef;
	$self->{"dna"}         = undef;
	my %temp_hash = ();
	$self->{"read_names"}        = \%temp_hash;
	$self->{"custom_chunk_size"} = 0;
	read_phylosift_config( self => $self );
	bless($self);
	return $self;
}

=head2 initialize
    
    Initializes the variables for the Phylosift object
    Using the standard pathnames and the filename

=cut

sub initialize {
	my $self        = shift;
	my %args        = @_;
	my $mode        = $args{mode};
	my $readsFile   = $args{file_1} || "";
	my $readsFile_2 = $args{file_2} || "";
	debug "READSFILE\t".$readsFile."\n" if length($readsFile_2);
	my $position = rindex( $readsFile, "/" );
	$self->{"fileName"}    = substr( $readsFile, $position + 1, length($readsFile) - $position - 1 );
	$self->{"workingDir"}  = getcwd;
	$self->{"mode"}        = $mode;
	$self->{"readsFile"}   = $readsFile;
	$self->{"readsFile_2"} = $readsFile_2;

	unless ( defined($Phylosift::Settings::file_dir) ) {
		$Phylosift::Settings::file_dir = $self->{"workingDir"}."/PS_temp/".$self->{"fileName"};
	}
	$self->{"blastDir"} = $Phylosift::Settings::file_dir."/blastDir";
	$self->{"alignDir"} = $Phylosift::Settings::file_dir."/alignDir";
	$self->{"treeDir"}  = $Phylosift::Settings::file_dir."/treeDir";
	$self->{"dna"}      = 0;
	%{ $self->{"read_names"} } = ();
	$self->{"run_info"} = Phylosift::Utilities::load_run_info( self => $self ) unless $self->{"mode"} eq 'all' && !$Phylosift::Settings::continue;

	# process defaults again now that command-line params have been parsed,
	# just in case any defaults depend on command-line settings
	Phylosift::FastSearch::set_default_values( post => 1 );
	Phylosift::MarkerAlign::set_default_values( post => 1 );
	Phylosift::pplacer::set_default_values( post => 1 );
	Phylosift::Summarize::set_default_values( post => 1 );

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

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Phylosift::Phylosift;

    my $foo = Phylosift::Phylosift->new();
    ...

=head1 SUBROUTINES/METHODS

=head2 run
    
    Runs the PhyloSift pipeline according to the functions passed as arguments

=cut

my $continue = 0;
my ( $mode, $readsFile, $readsFile_2, $fileName, $fileDir, $blastDir, $alignDir, $treeDir ) = "";
my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = 0;
my $workingDir = getcwd;

#where everything will be written when PhyloSift is running
sub run {
	my $self     = shift;
	my %args     = @_;
	my $force    = $args{force} || 0;
	my $custom   = $args{custom};
	my $continue = $args{cont}
	  || 0;    #continue is a reserved word, using a shortened version
	debug "force : $force\n";
	Phylosift::Utilities::print_citations();
	start_timer( name => "START" );

	read_phylosift_config( self => $self );
	run_program_check( self => $self );
	Phylosift::Utilities::data_checks( self => $self );

	#get shared lock for $marker_dir so nothing overrides the data while PS is running.
	my $lock_shared_markers = File::NFSLock->new( $Phylosift::Settings::marker_dir, LOCK_SH );
	my $lock_shared_ncbi    = File::NFSLock->new( $Phylosift::Settings::ncbi_dir,   LOCK_SH );
	file_check( self => $self );
	directory_prep( self => $self, force => $force, cont => $continue );

	# Forcing usage of updated markers
	debug "Using updated markers\n" if $Phylosift::Settings::updated;
	my @markers = Phylosift::Utilities::gather_markers( self => $self, marker_file => $custom );
	if ($Phylosift::Settings::extended) {
		@markers = Phylosift::Utilities::gather_markers( self => $self, path => $Phylosift::Settings::markers_extended_dir );
	}
	my $RUNINFO = ps_open( ">>".Phylosift::Utilities::get_run_info_file( self => $self ) );
	Phylosift::Summarize::print_run_info( self => $self, OUTPUT => $RUNINFO );
	close($RUNINFO);

	#changed > to >> so that we append to that file every time a run is restarted and we can keep track if there was a change in DB or in command
	#Also it prevents the file from being clobbered every time a run is restarted.
	if ( $self->{"mode"} eq 'search' || $self->{"mode"} eq 'all' ) {
		$self = run_search( self => $self, cont => $continue, marker => \@markers );
		$self->{"mode"} = 'align' if $continue;
	}
	debug "MODE :: ".$self->{"mode"}."\n";

	if ( defined($Phylosift::Settings::start_chunk) && defined($Phylosift::Settings::chunks) ) {
		for ( my $c = $Phylosift::Settings::start_chunk; $c < $Phylosift::Settings::start_chunk + $Phylosift::Settings::chunks; $c++ ) {
			run_later_stages( self => $self, cont => $continue, marker => \@markers, chunk => $c );
		}
	} else {
		run_later_stages( self => $self, cont => $continue, marker => \@markers );
	}
	$lock_shared_markers->unlock();
	$lock_shared_ncbi->unlock();
	Phylosift::Utilities::end_timer( name => "START" );
}

sub run_later_stages {
	my %args     = @_;
	my $self     = $args{self} || miss("self");
	my $continue = $args{cont};
	if ( $self->{"mode"} eq 'align' ) {
		$self = run_marker_align(@_);
		$self->{"mode"} = 'placer' if $continue;
	}
	if ( $self->{"mode"} eq 'placer' ) {
		$self = run_pplacer(@_);
		$self->{"mode"} = 'summarize' if $continue;
	}
	if ( $self->{"mode"} eq 'summarize' ) {
		$self = taxonomy_assignments(@_);
	}
}

=head2 read_phylosift_config

    Reads the Phylosift configuration file and assigns the file paths to the required directories

=cut

sub read_phylosift_config {
	my %args          = @_;
	my $custom_config = $Phylosift::Settings::configuration;
	debug "Reading config file : $custom_config\n" if $custom_config;

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
	if ( defined $Phylosift::Settings::CHUNK_MAX_SEQS ) {
		$args{self}->{"custom_chunk_size"} = 1;
	}

	#apply the command line parameters to override the RC files
}

=head2 fileCheck

    Checks if the files passed to the Phylosift object exist and are not empty

=cut

sub file_check {
	my %args = @_;
	my $self = $args{self} || miss("self");
	return $self if $self->{stdin};
	if ( !-e $self->{"readsFile"} ) {
		croak( $self->{"readsFile"}."  was not found \n" );
	}

	#check if the input file is a file and not a directory
	if ( !-f $self->{"readsFile"} ) {
		croak( $self->{"readsFile"}." is not a plain file, could be a directory\n" );
	}
	if ( $self->{"readsFile_2"} ne "" ) {

		#check the input file exists
		if ( !-e $self->{"readsFile_2"} ) {
			croak( $self->{"readsFile_2"}." was not found\n" );
		}

		#check if the input file is a file and not a directory
		if ( !-e $self->{"readsFile_2"} ) {
			croak( $self->{"readsFile_2"}." is not a plain file, could be a directory\n" );
		}
	}
	return $self;
}

=head2 run_program_check

    Runs a check on the programs that will be used through the pipeline to make sure they are
    available to the user and are the versions Phylosift was tested with.

=cut

sub run_program_check {
	my %args = @_;
	my $self = $args{self} || miss("self");

	#check if the various programs used in this pipeline are installed on the machine
	my $progCheck = Phylosift::Utilities::program_checks($self);
	if ( $progCheck != 0 ) {
		croak "A required program was not found during the checks aborting\n";
	} elsif ( $progCheck == 0 ) {
		debug "All systems are good to go, continuing the screening\n";
	}
	return $self;
}

=head2 directoryPrep

    Prepares the temporary Phylosift directory by deleting old runs and/or creating the correct directory structure
    
=cut

sub directory_prep {
	my %args     = @_;
	my $self     = $args{self} || miss("self");
	my $force    = $args{force};
	my $continue = $args{cont};
	##remove the directory from a previous run
	if ( -e $Phylosift::Settings::file_dir && $self->{"mode"} eq 'all' ) {
		Phylosift::Utilities::previous_run_handler( self => $self, dir => $Phylosift::Settings::file_dir, force => $force ) unless $continue;
		debug "Found existing directory, continuing a previous run\n" if $continue;
	}

	#create a directory for the Reads file being processed.
	`mkdir -p "$Phylosift::Settings::file_dir"`
	  unless ( -e $Phylosift::Settings::file_dir );
	debug "MODE : ".$self->{"mode"}."\n";
	`mkdir -p "$self->{"blastDir"}"` if ( !-e $self->{"blastDir"} && $self->{"mode"} =~ /search|all/ );
	`mkdir -p "$self->{"alignDir"}"` if ( !-e $self->{"alignDir"} && $self->{"mode"} =~ /align|all/ );
	`mkdir -p "$self->{"treeDir"}"`  if ( !-e $self->{"treeDir"}  && $self->{"mode"} =~ /place|all/ );
	return $self;
}

=head2 taxonomy_assignments

    performs the appropriate checks before running the taxonomy classification parts of the pipeline

=cut

sub taxonomy_assignments {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $continue    = $args{cont};
	my $markListRef = $args{marker} || miss("marker");
	my $chunk       = $args{chunk};
	Phylosift::Utilities::start_timer( name => "taxonomy assignments" );
	Phylosift::Summarize::summarize( self => $self, marker_reference => $markListRef, chunk => $chunk );
	Phylosift::Utilities::end_timer( name => "taxonomy assignments" );
	return $self;
}

=head2 compare

=cut

sub compare {
	my %args = @_;
	my $self = $args{self} || miss("self");
	debug "RUNNING Compare\n";
	Phylosift::Comparison::compare( self => $self, parent_dir => "./" );
}

=head2 run_pplacer

    Performs the appropriate checks before runing the Read placement parts of the pipeline

    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub run_pplacer {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $continue    = $args{cont} || 0;
	my $markListRef = $args{marker} || miss("marker");
	my $chunk       = $args{chunk};
	Phylosift::Utilities::start_timer( name => "runPPlacer" );
	Phylosift::pplacer::pplacer( self => $self, marker_reference => $markListRef, chunk => $chunk );
	Phylosift::Utilities::end_timer( name => "runPPlacer" );
	return $self;
}

=head2 run_marker_align

    Run the Hmm hit verification and the hit alignments to prep the hits for Pplacer
    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub run_marker_align {
	my %args     = @_;
	my $self     = $args{self} || miss("self");
	my $continue = $args{cont} || 0;
	my $markRef  = $args{marker} || miss("marker");
	my $chunk    = $args{chunk};
	Phylosift::Utilities::start_timer( name => "Alignments" );

	#clearing the alignment directory if needed
	my $alignDir = $self->{"alignDir"};

	#`rm "$alignDir"/*` if (<"$alignDir"/*>) && $Phylosift::Settings::force;

	#Align Markers
	my $threadNum = 1;
	Phylosift::MarkerAlign::MarkerAlign( self => $self, marker_reference => $markRef, chunk => $chunk );
	Phylosift::Utilities::end_timer( name => "Alignments" );
	return $self;
}

=head2 run_search

    Searches reads against the various databases of the pipeline
    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub run_search {
	my %args          = @_;
	my $self          = $args{self} || miss("self");
	my $continue      = $args{cont} || 0;
	my $markerListRef = $args{marker} || miss("marker");
	Phylosift::Utilities::start_timer( name => "runBlast" );

	#clearing the blast directory
	#my $blastDir = $self->{"blastDir"};
	#`rm "$blastDir"/*` if (<"$blastDir"/*>);

	#run Searches
	Phylosift::FastSearch::run_search( self => $self, marker_reference => $markerListRef );
	Phylosift::Utilities::end_timer( name => "runBlast" );
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

=over

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
