package Phylosift::Utilities;
use strict;
use warnings;
use FindBin qw($Bin);
BEGIN { unshift( @INC, "$FindBin::Bin/../legacy/" ) if $] < 5.01; }
use File::Basename;
use File::NFSLock qw(uncache);
use Fcntl qw(LOCK_EX);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Phylosift::Settings;
use Bio::Align::Utilities qw(:all);
use Bio::TreeIO;
use Bio::Tree::Tree;
use POSIX ();
use Carp;
use Cwd;

our $VERSION = "v1.0.1";

use Exporter;
use vars qw[ @EXPORT @EXPORT_OK %EXPORT_TAGS @ISA ];
@ISA       = 'Exporter';
@EXPORT    = qw[start_timer end_timer debug miss ps_open];
@EXPORT_OK = qw[];
%EXPORT_TAGS = (
				 STD => \@EXPORT,
				 all => [ @EXPORT, @EXPORT_OK ],
);
our $debuglevel = 0;
my %timers;

sub debug {
	my $msg = shift;
	my $msglevel = shift || 1;
	print $msg if $debuglevel >= $msglevel;
}

=head1 NAME

Phylosift::Utilities - Implements miscellaneous accessory functions for Phylosift

=head1 VERSION

Version 0.01

=cut

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
=cut

sub miss {
	my $arg_name = shift;
	$Carp::Verbose = 1;
	croak("Missing required argument $arg_name");
}

sub ps_open {
	my $open_string = shift;
	my $result = open( my $FH, $open_string );
	unless ($result) {
		my $type = "read from";
		$type = "write to"  if $open_string =~ /^>/;
		$type = "append to" if $open_string =~ /^>>/;
		$type = "append to" if $open_string =~ /^\+>/;
		$type = "pipe"      if $open_string =~ /\|$/;
		$Carp::Verbose = 1;
		croak("Unable to $type $open_string");
	}
	return $FH;
}

sub get_program_path {
	my %args      = @_;
	my $progname  = $args{prog_name};
	my $progpath  = $args{prog_path};
	my $progcheck = "";

	# if a path was already given try that first
	if (    defined($progpath)
		 && $progpath ne ""
		 && -x $progpath."/".$progname )
	{
		$progcheck = $progpath."/".$progname;
	}

	# then check the directories from where the script is running
	$progcheck = $Bin."/".$progname
	  unless ( $progcheck =~ /$progname/ || !( -x $Bin."/".$progname ) );
	$progcheck = $Bin."/bin/".$progname
	  unless ( $progcheck =~ /$progname/
			   || !( -x $Bin."/bin/".$progname ) );

	# then check the OS and use Mac binaries if needed
	if ( $^O =~ /arwin/ ) {
		$progcheck = $Bin."/../osx/".$progname
		  unless ( $progcheck =~ /$progname/
				   && !( -x $Bin."/".$progname ) );
		$progcheck = $Bin."/../osx/".$progname
		  if ( $progcheck =~ /$Bin\/bin/ );    # don't use the linux binary!
	}

	# last ditch attempt, try the system path and hope we get the right version
	unless ( -x $progcheck ) {
		$progcheck = `which "$progname" 2> /dev/null`;
		chomp $progcheck;
	}

	return $progcheck;
}

# external programs used by Phylosift

our %marker_lookup = ();

=head2 program_checks

checks the program requirements for PhyloSift
writes to STDERR if a program is missing or the wrong version is installed
returns 1 or 0 depending on success of failure.

=cut

sub program_checks {
	eval 'require Bio::Seq;';
	if ($@) {
		carp "Bioperl was NOT found\n";
		return 1;
	}
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::pplacer,
									  value     => get_program_path(
																 prog_name => "pplacer",
																 prog_path => $Phylosift::Settings::pplacer_path
									  )
	);
	if ( $Phylosift::Settings::pplacer eq "" ) {

		#program not found return;
		carp("pplacer not found");
		return 1;
	} else {
		`$Phylosift::Settings::pplacer --version` =~ m/v1.1.alpha(\d+)/;
		if ( defined($1) && $1 < 14 ) {

			# pplacer was found but the version doesn't match the one tested with Phylosift
			carp("Warning : a different version of pplacer was found. PhyloSift was tested with pplacer v1.1.alpha14\n");
		}
	}
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::guppy,
									  value     => get_program_path(
																 prog_name => "guppy",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::taxit,
									  value     => get_program_path(
																 prog_name => "taxit",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::rppr,
									  value     => get_program_path(
																 prog_name => "rppr",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::cmalign,
									  value     => get_program_path(
																 prog_name => "cmalign",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::hmmalign,
									  value     => get_program_path(
																 prog_name => "hmmalign",
																 prog_path => $Phylosift::Settings::hmmer3_path
									  )
	);
	if ( $Phylosift::Settings::hmmalign eq "" ) {

		#program not found return;
		carp("HMMER3 not found");
		return 1;
	} elsif ( `$Phylosift::Settings::hmmalign -h` !~ m/HMMER 3.0/ ) {

		# pplacer was found but the version doens't match the one tested with Phylosift
		carp "Warning : a different version of HMMER was found. PhyloSift was tested with HMMER 3.0rc1\n";
	}
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::hmmbuild,
									  value     => get_program_path(
																 prog_name => "hmmbuild",
																 prog_path => $Phylosift::Settings::hmmer3_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::readconciler,
									  value     => get_program_path(
																 prog_name => "readconciler",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::segment_tree,
									  value     => get_program_path(
																 prog_name => "segment_tree",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::pda,
									  value     => get_program_path(
																 prog_name => "pda",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::fasttree,
									  value     => get_program_path(
																 prog_name => "FastTree",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::lastdb,
									  value     => get_program_path(
																 prog_name => "lastdb",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::lastal,
									  value     => get_program_path(
																 prog_name => "lastal",
																 prog_path => $Phylosift::Settings::ps_path
									  )
	);

	return 0;
}

=head2 read_marker_summary

reads the marker_summary file from a specific directory

=cut

sub read_marker_summary {
	my %args        = @_;
	my $self        = $args{self} || miss("PS object");
	my $directory   = $args{path} || miss("marker_summary Path");
	my %return_hash = ();
	return \%return_hash unless -e "$directory/marker_summary.txt";
	my $IN     = ps_open( $directory."/marker_summary.txt" );
	my $header = <$IN>;                                         #read header line
	while (<$IN>) {
		chomp($_);
		my @line = split( /\t/, $_ );
		next if $line[0] eq 'total';
		$return_hash{ $line[0] } = $line[1];
	}
	return \%return_hash;
}

=head2 print_bat_signal

Prints the batsignal. Should only be used in case of emergency

=cut

sub print_bat_signal {
	my $signal = "       _==/          i     i          \\==_\n";
	$signal .= "     /XX/            |\\___/|            \\XX\n";
	$signal .= "   /XXXX\\            |XXXXX|            /XXXX\\\n";
	$signal .= "  |XXXXXX\\_         _XXXXXXX_         _/XXXXXX|\n";
	$signal .= " XXXXXXXXXXXxxxxxxxXXXXXXXXXXXxxxxxxxXXXXXXXXXXX\n";
	$signal .= "|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|\n";
	$signal .= "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
	$signal .= "|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX|\n";
	$signal .= " XXXXXX/^^^^^\\XXXXXXXXXXXXXXXXXXXXX/^^^^^\\XXXXXX\n";
	$signal .= "  |XXX|       \\XXX/^^\\XXXXX/^^\\XXX/       |XXX|\n";
	$signal .= "    \\XX\       \\X/    \\XXX/    \\X/       /XX/\n";
	$signal .= "       \"\\       \"      \\X/      \"      /\"\n";
	print STDERR $signal;
}

=head2 print_marker_summary

reads the marker_summary file from a specific directory

=cut

sub print_marker_summary {
	my %args        = @_;
	my $self        = $args{self} || miss("PS object");
	my $dir         = $args{path} || miss("marker_summary path");
	my $summary_ref = $args{summary} || miss("Summary hash");
	my $OUT         = ps_open( ">".$dir."/marker_summary.txt" );
	my %summary     = %{$summary_ref};
	my $total       = 0;
	foreach my $marker ( sort keys %summary ) {
		print $OUT $marker."\t".$summary{$marker}."\n" if $summary{$marker} > 0;
		$total += $summary{$marker} unless $marker eq 'concat';
	}
	print $OUT "total\t$total\n";
	close($OUT);
}

=head2 dataChecks

Check for requisite PhyloSift marker datasets

=cut

sub get_data_path {
	my %args      = @_;
	my $dataname  = $args{data_name};
	my $datapath  = $args{data_path};
	my $datacheck = "";
	if ( defined($datapath) && $datapath ne "" ) {
		$datacheck = $datapath."/".$dataname;
	} else {
		my $scriptpath = dirname($0);
		$scriptpath =~ s/bin\/?$//g;
		$datacheck = $scriptpath."/share/phylosift/".$dataname;
		return $datacheck if ( -x $datacheck );

		# if the system data dir doesn't exist, default to the user's home
		$datacheck = $ENV{"HOME"}."/share/phylosift/".$dataname;
	}
	return $datacheck;
}

sub get_marker_version {
	my %args = @_;
	my $path = $args{path} || miss("path");
	return unless -e "$path/version.txt";
	my $MFILE = ps_open("$path/version.txt");
	my $url   = <$MFILE>;
	chomp $url;
	my $timestamp = <$MFILE>;
	chomp $timestamp;
	return ( $url, $timestamp );
}

=head2 download_data

Downloads the data given a url and a destination for the data

=cut

sub download_data {
	my %args        = @_;
	my $url         = $args{url};
	my $destination = $args{destination};

	# FIXME this is insecure!
	# but then again, so is just about every other line of code in this program...
	`mkdir -p "$destination"; cd $destination/.. ; curl -LO $url`;
	debug "URL : $url\n";
	$url =~ /\/(\w+)\.tgz/;
	my $archive = $1;
	debug "ARCHIVE : $archive\n";
	if ( -e "$destination/.. " ) {
		`rm -rf "$destination/.."`;
	}
	`cd "$destination/../" ; tar xzf $archive.tgz ; touch $archive`;
	`rm "$destination/../$archive.tgz"`;
}

Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::marker_base_url,
								  value     => "http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift_markers" );
Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::ncbi_url,
								  value     => "http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz" );

=head2 marker_update_check

Checks with the remote location to see if the markers have been updated since the last download
If the downloadable marker package is newer than what is currently being used then download the newer version.
=cut

sub marker_update_check {
	my %args        = @_;
	my $self        = $args{self};
	my $url         = $args{url};
	my $marker_path = $args{dir};
	my $skip_lock   = $args{skip_lock} || 0;
	debug "Skipping check for marker updates on server\n"
	  if $Phylosift::Settings::disable_update_check && !$self->{"removed_markers"};
	return
	  if $Phylosift::Settings::disable_update_check && !$self->{"removed_markers"};    # bail out if we're not supposed to be here
	                                                                                   #return
	                                                                                   #  if defined($Phylosift::Settings::marker_update_check)
	                                                                                   #  && $Phylosift::Settings::marker_update_check == 0;

	eval {
		require LWP::Simple;
		LWP::Simple->import();
	};
	my $result          = $@;
	my $get_new_markers = 0;
	my $modified_time   = 0;
	$url =~ s/\/\//\//g;
	$url =~ s/http:/http:\//g;

	unless ($result) {

		my ( $content_type, $document_length, $modified, $expires, $server ) = head($url);
		$modified_time = $modified;
		debug "MARKER_PATH : ".$marker_path."\nURL : $url\n" if $skip_lock;

		if ( -x $marker_path ) {
			my ( $m_url, $m_timestamp ) = get_marker_version( path => $marker_path );
			if ( defined($m_url) ) {
				$m_url =~ s/\/\//\//g;
				$m_url =~ s/http:/http:\//g;
				debug "TEST LOCAL :".localtime($m_timestamp)."\n" if $skip_lock;
				if ( !defined($modified_time) ) {
					warn "Warning: unable to connect to marker update server, please check your internet connection\n" if $skip_lock;
				} elsif ( $modified_time > $m_timestamp ) {
					debug "TEST REMOTE:".localtime($modified_time)."\n" if $skip_lock;
					#only warning to be displayed in the first iteration of marker_update_check since we don't want to come back through it.
					warn
					  "A newer version of the marker data exists on the server. Move or remove the current marker DB at $marker_path to trigger a download.\n"
					  unless $skip_lock;
					$get_new_markers = 0; # we don't want to flag to get new markers, let the user remove their DB manually.
				} elsif ( $url ne $m_url ) {
					warn "The marker update URL differs from the local marker DB copy, updating" if $skip_lock;
					warn "local url $m_url\n"                                                    if $skip_lock;
					warn "update url $url\n"                                                     if $skip_lock;
					$get_new_markers = 1;
				}
			} else {
				warn "Marker version unknown, downloading markers\n" if $skip_lock;
				$get_new_markers = 1;
			}
		} else {
			if ( !defined($modified_time) ) {
				croak
				  "Marker data not found and unable to connect to marker update server, please check your phylosift configuration and internet connection!\n"
				  if $skip_lock;
			}
			warn "Unable to find marker data!\n" if $skip_lock;
			$get_new_markers = 1;
		}
	} else {
		warn "Unable to check for updates because perl module LWP::Simple is not available. Please use cpan to install it.\n" if $skip_lock;
		$get_new_markers = 1 unless -e "$marker_path/version.txt";
	}
	my $lock_excl_markers;
	if ($get_new_markers) {
		`mkdir -p $marker_path`;
		warn "Attempting to get exclusive rights to the marker DB. Possibly waiting for another PhyloSift to release the file lock\n" unless $skip_lock;
		# Getting exclusive rights to the marker DB.
		# Rights are released when we return from the function (where lock_excl_markers is declared)
		# Skip_lock is used to let the function know to only print some error messages once. When the DB is tested when the lock is exclusive.
		# This implies that when starting many PS at once without a DB will most lilely result in all PS running sequentially.
		$lock_excl_markers = File::NFSLock->new( $marker_path, LOCK_EX ) unless $skip_lock;
		$get_new_markers = marker_update_check(
												self      => $self,
												dir       => $marker_path,
												url       => $url,
												skip_lock => 1
		) unless $skip_lock;
		warn "Downloading from $url\n" if $skip_lock;
		download_data( url => $url, destination => $marker_path ) if $skip_lock;
		return $get_new_markers if $skip_lock;
		my $VOUT = ps_open(">$marker_path/version.txt");
		print $VOUT "$url\n";
		print $VOUT "$modified_time\n";
	}
	my $indexed_markers = 0;
	if ( $get_new_markers || $self->{"removed_markers"} ) {
		my @markers = gather_markers( self => $self, path => $marker_path, force_gather => 1 );
		$indexed_markers = index_marker_db( self => $self, markers => \@markers, path => $marker_path );
	}
	return $indexed_markers;
}

=head2 data_checks

Performs checks to see if the NCBI data available for download is newer than what is currently being used by PhyloSift
If the available data is newer than the used data, download the newer version.

=cut

sub data_checks {
	my %args = @_;
	my $self = $args{self};

	#
	# determine where to look for markers
	my $mbase;
	if ( defined($Phylosift::Settings::marker_base_url) ) {
		$mbase = $Phylosift::Settings::marker_base_url;
	}
	if ( defined($Phylosift::Settings::marker_url) ) {
		$mbase = $Phylosift::Settings::marker_url;
	}
	my $marker_update_url           = "$mbase/markers.tgz";
	my $markers_extended_update_url = "$mbase/markers_extended.tgz";
	my $indexed_markers             = 0;

	#
	# first check for the standard marker directory
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::marker_dir,
									  value     => get_data_path(
															  data_name => "markers",
															  data_path => $Phylosift::Settings::marker_path
									  )
	);

	#checking if the PMMPROK markers are present at the same time as DNGNGWU markers
	my @pmprok  = glob("$Phylosift::Settings::marker_dir/PMPROK*");
	my @dngngwu = glob("$Phylosift::Settings::marker_dir/DNGNGWU*");
	if ( scalar(@pmprok) > 0 && scalar(@dngngwu) > 0 ) {
		warn("Old AND new core markers detected. Removing the old version\nForcing a re-index of the reference database...\n");
		my $cmd = "rm -rf $Phylosift::Settings::marker_dir/PMPROK*";
		`$cmd`;
		$self->{"removed_markers"} = 1;
	}
	$indexed_markers = marker_update_check(
											self => $self,
											dir  => $Phylosift::Settings::marker_dir,
											url  => $marker_update_url
	);

	#
	# now look for the extended marker directory
	Phylosift::Settings::set_default(
									  parameter => \$Phylosift::Settings::markers_extended_dir,
									  value     => get_data_path(
															  data_name => "markers_extended",
															  data_path => $Phylosift::Settings::extended_marker_path
									  )
	);

	if ($Phylosift::Settings::extended) {
		$indexed_markers = marker_update_check(
												self => $self,
												dir  => $Phylosift::Settings::markers_extended_dir,
												url  => $markers_extended_update_url
		);
	}

	#
	# now check for the NCBI taxonomy data
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::ncbi_dir,
									  value     => get_data_path( data_name => "ncbi", data_path => $Phylosift::Settings::ncbi_path ) );
	debug "Skipping check for NCBI updates on server\n" if $Phylosift::Settings::disable_update_check;
	return $indexed_markers if $Phylosift::Settings::disable_update_check;
	ncbi_update_check(
					   self => $self,
					   dir  => $Phylosift::Settings::ncbi_dir,
					   url  => $Phylosift::Settings::ncbi_url
	);
	return $indexed_markers;
}

=head2 ncbi_update_check

checks to see if a new version of the ncbi data is available on the server.
Download it if needed.

=cut

sub ncbi_update_check {
	my %args     = @_;
	my $self     = $args{self};
	my $url      = $args{url};
	my $ncbi_dir = $args{dir};
	$ncbi_dir =~ m/^(\S+)\/[^\/]+$/;    #Need to work from the parent directory from ncbi_dir
	my $ncbi_path = $1;
	my $skip_lock = $args{skip_lock};
	my $have_lwp  = 0;
	my @lwp_result;
	eval {
		require LWP::Simple;
		LWP::Simple->import();
		$have_lwp   = 1;
		@lwp_result = LWP::Simple::head("$url");
	  }
	  or do {
	  };
	my ( $content_type, $document_length, $modified_time, $expires, $server ) = @lwp_result;

	if ( -x $ncbi_dir ) {
		my $ncbi_time = ( stat($ncbi_dir) )[9];
		if ( !defined($modified_time) ) {
			warn "Warning: unable to connect to NCBI taxonomy update server, please check your internet connection\n" if $skip_lock && $have_lwp;
		} elsif ( $modified_time > $ncbi_time ) {
			warn "Found newer version of NCBI taxonomy data!\n" if $skip_lock;
			warn "Downloading from $url\n"                      if $skip_lock;
			my $lock_ex;
			$lock_ex = File::NFSLock->new( $ncbi_path, LOCK_EX ) unless $skip_lock;
			ncbi_update_check( self => $self, dir => $ncbi_dir, url => $url, skip_lock => 1 ) unless ($skip_lock);
			download_data( url => $url, destination => $ncbi_dir ) if $skip_lock;
			$lock_ex->unlock() unless $skip_lock;
		}
	} else {
		if ( !defined($modified_time) && $have_lwp ) {
			croak "NCBI taxonomy data not found and unable to connect to update server, please check your phylosift configuration and internet connection!\n"
			  if $skip_lock;
		}
		warn "Unable to find NCBI taxonomy data!\n" if $skip_lock;
		warn "Downloading from $url\n"              if $skip_lock;
		my $lock_ex;
		$lock_ex = File::NFSLock->new( $ncbi_path, LOCK_EX ) unless $skip_lock;
		ncbi_update_check( self => $self, dir => $ncbi_dir, url => $url, skip_lock => 1 ) unless ($skip_lock);
		download_data( url => $url, destination => $ncbi_dir ) if $skip_lock;
		return if $skip_lock;
		$lock_ex->unlock() unless $skip_lock;
	}
}

=head2 previous_run_handler

Exits and Prints out the error statement when a previous run was detected and -f was not specified
OR
Removes the appropriate directory when -f was specified

=cut

sub previous_run_handler {
	my %args  = @_;
	my $self  = $args{self} || miss("PS Ojbect in Utilities::previous_run_handler");
	my $force = $args{force};
	my $dir   = $args{dir} || miss("dir in Utilities::previous_run_handler");
	if ($force) {
		debug( "deleting an old run\n", 0 );
		debug "$dir\n";
		`rm -rf "$dir"`;
	} else {
		croak(  "A previous run was found using the same output name. Aborting the current run.\n"
			   ."Either delete that run from "
			   .$dir
			   .", or force overwrite with the -f command-line option\n" );
	}
}

=head2 fasta2stockholm

Convert a bunch of fasta files to stockholm format
This code is adapted from fasta2stockholm.pl, (c) Ian Holmes and licensed under the GPL
See https://github.com/ihh/dart/ for the original source code

=cut

sub fasta2stockholm {
	my %args     = @_;
	my $fasta    = $args{fasta};
	my $output   = $args{output};
	my $STOCKOUT = ps_open(">$output");

	# read FASTA file
	my @seq;
	my @name;
	my $name;
	my $curseq = "";
	my $FASTA  = ps_open("<$fasta");
	while (<$FASTA>) {
		if (/^\s*>\s*(\S+)/) {
			if ( length($curseq) > 0 ) {
				push @seq,  $curseq;
				push @name, $name;
				$curseq = "";
			}
			$name = $1;
		} else {
			if ( /\S/ && !defined $name ) {
				warn "Ignoring: $_";
			} else {
				s/\s//g;
				$curseq .= $_;
			}
		}
	}
	if ( length($curseq) > 0 ) {
		push @seq,  $curseq;
		push @name, $name;
		$curseq = "";
	}
	close $FASTA;

	# check all seqs are same length
	my $length;
	my $lname;
	for ( my $sI = 0; $sI < @name; $sI++ ) {
		my $nname = $name[$sI];
		my $sseq  = $seq[$sI];
		my $l     = length $sseq;
		if ( defined $length ) {
			croak "Sequences not all same length ($lname is $length, $nname is $l)"
			  unless $length == $l;
		} else {
			$length = length $sseq;
			$lname  = $nname;
		}
	}

	# print Stockholm output
	print $STOCKOUT "# STOCKHOLM 1.0\n";
	for ( my $sI = 0; $sI < @name; $sI++ ) {
		print $STOCKOUT $name[$sI], " ", $seq[$sI], "\n";
	}
	print $STOCKOUT "//\n";
}

=head2 stockholm2fasta

Convert a stockholm file to fasta format
This code is adapted from stockholm2fasta.pl, (c) Ian Holmes and licensed under the GPL
See https://github.com/ihh/dart/ for the original source code

=cut

sub stockholm2fasta {
	my %args     = @_;
	my $columns  = $args{columns} || 80;    # number of columns in fasta
	my $gapped   = $args{gapped} || 1;      # should gapped fasta (aligned) be written?
	my $sorted   = $args{sorted} || 1;      # should sequences be sorted?
	my $STREAMIN = $args{in};
	my @seq;
	my @names;
	my $outbuffer = "";
	my $seq_line  = 0;                      # counter for which line we're on in each block.

	while ( my $line = <$STREAMIN> ) {
		if ( $line !~ /\S/ || $line =~ /^\s*#/ || length($line) < 2 ) {
			$seq_line = 0;
			next;
		}
		if ( $line =~ /^\s*\/\// ) {
			$outbuffer .= printseq(
									columns => $columns,
									sorted  => $sorted,
									seq     => \@seq,
									names   => \@names,
									out     => $args{out}
			);
		} else {
			chomp $line;
			my ( $name, $seq ) = split /\s+/, $line;
			$seq =~ s/[\.\-]//g unless $gapped;
			$seq[$seq_line] .= $seq;
			$names[$seq_line] = $name;
			$seq_line++;
		}
	}
	return $outbuffer;
}

sub printseq {
	my %args  = @_;
	my @seq   = @{ $args{seq} };
	my @names = @{ $args{names} };
	my $out   = "";
	for ( my $j = 0; $j < @names; $j++ ) {
		$out .= ">$names[$j]\n";
		for ( my $i = 0; $i < length( $seq[$j] ); $i += $args{columns} ) {
			$out .= substr( $seq[$j], $i, $args{columns} )."\n";
		}
	}
	return $out;
}

=head2 read_name_table

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub read_name_table {
	my %args      = @_;
	my $markerDir = $args{marker_directory};
	my %result;
	my $ALINAMES = ps_open("grep \">\" $markerDir/*.ali |");
	while ( my $line = <$ALINAMES> ) {
		$line =~ s/.+\:\>//g;
		my $commonName;
		if ( $line =~ /\{(.+)\}.+\[.+\]/ ) {
			$commonName = $1;
		} elsif ( $line =~ /\[(.+)\]$/ ) {
			$commonName = $1;
		}
		my @fields = split( /\s+/, $line );
		$result{ $fields[0] } = $commonName;
	}
	return %result;
}

=head2 get_marker_length

reads in a marker name and returns the number of sequences used to build the HMM profile or CM

=cut

sub get_marker_length {
	my %args   = @_;
	my $self   = $args{self};
	my $marker = $args{marker};
	my $length = 0;
	if ( is_protein_marker( marker => $marker ) ) {
		my $hmm_file = get_marker_hmm_file( self => $self, marker => $marker );
		return unless -e $hmm_file;
		my $HMM = ps_open($hmm_file);
		while ( my $line = <$HMM> ) {
			if ( $line =~ /LENG\s+(\d+)/ ) {
				$length = $1;
				last;
			}
		}
	} else {
		my $cm_file = get_marker_cm_file( self => $self, marker => $marker );
		my $CM = ps_open($cm_file);
		while ( my $line = <$CM> ) {
			if ( $line =~ /CLEN\s+(\d+)/ ) {
				$length = $1;
				last;
			}
		}
	}
	return $length;
}

=head2 make_dummy_file

Creates a dummy file

=cut

sub make_dummy_file {
	my %args          = @_;
	my $self          = $args{self};
	my $marker        = $args{marker};
	my $gapmultiplier = $args{gap_multiplier};
	my $len           = get_marker_length( self => $self, marker => $marker );
	my $glen;
	$glen = "P" x $len x $gapmultiplier if $gapmultiplier == 1;
	$glen = "A" x $len x $gapmultiplier if $gapmultiplier == 3;
	my $newseq = Bio::LocatableSeq->new(
										 -seq      => $glen,
										 -alphabet => 'dna',
										 -id       => "dummydummydummy",
										 start     => 1,
										 end       => ( $len * $gapmultiplier )
	);
	my $aln = Bio::SimpleAlign->new();
	$aln->add_seq($newseq);
	return $aln;
}

=head2 get_marker_path

Determines the filesystem path to a marker. Searches for the marker in the base marker directory, the extended markers, and any local markers

=cut

sub get_marker_path {
	my %args = @_;
	my $marker = $args{marker} || miss("marker");

	# check for old-style marker first
	return "$Phylosift::Settings::marker_dir" if ( -e "$Phylosift::Settings::marker_dir/$marker.faa" );

	# check for new-style in standard directory
	return "$Phylosift::Settings::marker_dir" if ( -d "$Phylosift::Settings::marker_dir/$marker" );

	# check for new-style in extended directory
	return "$Phylosift::Settings::markers_extended_dir" if ( -d "$Phylosift::Settings::markers_extended_dir/$marker" );

	# TODO: check any local marker repositories
	warn "Could not find repository for marker $marker\n";

	#	$Carp::Verbose = 1;
	#	croak("Could not find repository for marker $marker\n");
}

=head2 get_marker_basename

Returns the base name of the marker -- the marker name without any directories prepended

=cut

sub get_marker_basename {
	my %args   = @_;
	my $marker = $args{marker};
	$marker =~ s/^.+\///g;
	return $marker;
}

=head2 get_marker_basename

Returns the full name of the marker -- the marker name with any directories prepended

=cut

sub get_marker_fullname {
	my %args   = @_;
	my $marker = $args{marker};
	return $marker_lookup{$marker};
}

=head2 get_alignment_marker_file 

Returns the alignment file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_alignment_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.ali";
	} else {
		return "$marker_path/$bname.updated.fasta";
	}
}

=head2 get_marker_aln_file

Returns the aligned fasta file for the marker

=cut

sub get_marker_aln_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.ali" if ( -e "$marker_path/$bname.ali" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.aln"
		  if ( -e "$marker_path/$marker/$bname.aln" );
		return "$marker_path/$marker/$bname.masked";
	} else {
		return "$marker_path/$marker.updated/$bname.aln"
		  if ( -e "$marker_path/$marker.updated/$bname.aln" );
		return "$marker_path/$marker.updated/$bname.masked";
	}
}

=head2 get_marker_rep_file

Returns the fasta file of unaligned full length representative sequences for the marker

=cut

sub get_marker_rep_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $updated     = $args{updated} || 0;
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $updated == 0 ) {
		return "$marker_path/$bname.faa" if ( -e "$marker_path/$bname.faa" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.rep";
	} else {
		return "$marker_path/$marker.updated/$bname.reps";
	}
}

=head2 get_marker_hmm_file

Returns the HMM file for the marker

=cut

sub get_marker_hmm_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $local       = $args{loc} || 0;                                       #using loc instead of local (reserved word)
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	return "$marker_path/$bname.hmm" if ( -e "$marker_path/$bname.hmm" );

	# using new-style marker directories
	return "$marker_path/$marker/$bname.hmm";
}

=head2 get_marker_cm_file

Returns the CM (infernal covarion model) file for the marker

=cut

sub get_marker_cm_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	return "$marker_path/$marker/$bname.cm";
}

=head2 get_marker_stockholm_file

Returns the stockholm file for the marker

=cut

sub get_marker_stockholm_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	return "$marker_path/$bname.stk" if ( -e "$marker_path/$bname.ali" );

	# using new-style marker directories
	return "$marker_path/$marker/$bname.stk";
}

=head2 get_marker_taxon_map

Returns the path to the lookup table between marker gene IDs and their taxa

=cut

sub get_marker_taxon_map {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	my $decorated   = get_decorated_marker_name(%args);
	my $deco        = "$marker_path/$decorated/$decorated.taxonmap";
	return $deco if -e $deco;
	return "$marker_path/$marker.updated/$bname.updated.taxonmap"
	  if $Phylosift::Settings::extended && -e "$marker_path/$marker.updated/$bname.updated.taxonmap";
	return "$marker_path/$marker/$bname.taxonmap" if $Phylosift::Settings::extended;
	return "$marker_path/$bname/$bname.taxonmap";
}

sub get_decorated_marker_name {
	my %args       = @_;
	my $name       = $args{marker} || miss("marker");
	my $dna        = $args{dna} || 0;
	my $updated    = $args{updated} || 0;
	my $sub_marker = $args{sub_marker};
	my $base       = $args{base} || 0;
	my $long       = $args{long} || 0;
	my $short      = $args{short} || 0;
	$name = get_marker_basename( marker => $name ) if $base;
	$name .= ".long"  if $long;
	$name .= ".short" if $short;
	$name .= ".codon" if $dna;
	##$name .= ".updated"
	##  if ( $updated || $Phylosift::Settings::updated ) && !$Phylosift::Settings::extended && $name !~ /^1[68]s/;
	$name .= ".updated"
	  if ( $updated || $Phylosift::Settings::updated ) && $name !~ /^1[68]s/;

	# rna markers don't get pruned
	$name .= ".sub".$sub_marker if defined($sub_marker);
	return $name;
}

=head2 get_gene_id_file

Returns the gene id file name

=cut

sub get_gene_id_file {
	my %args        = @_;
	my $marker      = $args{marker} || miss("marker");
	my $dna         = $args{dna} || 0;
	my $marker_path = get_marker_path( marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	my $decorated   = get_decorated_marker_name(%args);
	my $deco        = "$marker_path/$decorated/$decorated.gene_map";
	return $deco if -e $deco;
	return "$marker_path/$marker/$bname.gene_map" if $Phylosift::Settings::extended && -e "$marker_path/$marker/$bname.gene_map";
	return "$marker_path/$bname/$bname.gene_map" if -e "$marker_path/$bname/$bname.gene_map";
	return "$Phylosift::Settings::marker_dir/gene_ids.codon.txt" if $dna;
	return "$Phylosift::Settings::marker_dir/gene_ids.aa.txt";
}

=head2 is_protein_marker

Returns 1 if the marker named in 'marker' is protein sequence (0 if RNA or DNA)

=cut

sub is_protein_marker {
	my %args        = @_;
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );

	# only protein markers have HMMs
	return 0 if ( -e "$marker_path/$marker/$bname.cm" );
	return 1 if ( -e "$marker_path/$marker/$bname.hmm" );
	return 1 if ( -e "$marker_path/$bname.hmm" );
	return 1;
}

=head2 get_fasta_marker_file

Returns the fasta file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_fasta_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.faa";
	} else {
		return "$marker_path/$bname.faa";
	}
}

=head2 get_marker_package

Returns the path to the marker package

=cut

sub get_marker_package {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $decorated   = get_decorated_marker_name(%args);
	return "$marker_path/$decorated";
}

=head2 get_aligner_output_fasta
Returns the FastA file containing DNA read or contig alignments to the marker
given by markerName
=cut

sub get_aligner_output_fasta {
	my %args      = @_;
	my $marker    = $args{marker};
	my $chunk     = $args{chunk};
	my $chunky    = defined($chunk) ? ".$chunk" : "";
	my $decorated = get_decorated_marker_name( %args, base => 1 );
	return "$decorated$chunky.fasta";
}

=head2 get_read_placement_file
Returns the read placement Jplace file to the marker
given by markerName
=cut

sub get_read_placement_file {
	my %args      = @_;
	my $self      = $args{self};
	my $marker    = $args{marker};
	my $chunk     = $args{chunk};
	my $bname     = get_marker_basename( marker => $marker );
	my $decorated = get_decorated_marker_name( %args, base => 1 );
	return "$decorated.$chunk.jplace";
}

=head2 get_trimfinal_marker_file

Returns the .trimfinal file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_trimfinal_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.trimfinal"
		  if -e "$marker_path/$bname.trimfinal";
		return "$marker_path/$marker/$bname.masked"
		  if -e "$marker_path/$marker/$bname.masked";
		return "$marker_path/$marker/$bname.clean"
		  if -e "$marker_path/$marker/$bname.clean";
		return "$marker_path/$marker/$bname.aln";
	} else {
		return "$marker_path/$bname.trimfinal"
		  if -e "$marker_path/$bname.trimfinal";
		return "$marker_path/$marker/$bname.masked"
		  if -e "$marker_path/$marker/$bname.masked";
		return "$marker_path/$marker/$bname.clean"
		  if -e "$marker_path/$marker/$bname.clean";
		return "$marker_path/$marker/$bname.aln";
	}
}

=head2 get_trimfinal_fasta_marker_file

Returns the .trimfinal.fasta file for the markerName passed in as an argument
If the user chooses the updated markers, the updated file is returned instead

=cut

sub get_trimfinal_fasta_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.trimfinal.fasta"
		  if -e "$marker_path/$bname.trimfinal.fasta";
		return "$marker_path/$marker/$bname.aln";
	} else {
		return "$marker_path/$bname.trimfinal.fasta"
		  if -e "$marker_path/$bname.trimfinal.fasta";
		return "$marker_path/$marker/$bname.aln";
	}
}

=head2 get_tree_marker_file

Returns the .final.tre file from the marker directory 
The user chooses the updated or stock version

=cut

sub get_tree_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.final.tre"
		  if -e "$marker_path/$bname.final.tre";
		return "$marker_path/$marker/$bname.tre";
	} else {
		return "$marker_path/$bname.updated.tre";
	}
}

=head2 get_tree_stats_marker_file

Return the updated or stock version of the Tree stats file

=cut

sub get_tree_stats_marker_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.in_phyml_stats.txt";
	} else {
		return "$marker_path/$bname.updated.RAxML_info";
	}
}

=head2 get_ncbi_map_file

Returns the updated of stock version of the NCBI map file

=cut

sub get_ncbi_map_file {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $Phylosift::Settings::updated == 0 ) {
		return "$marker_path/$bname.ncbimap";
	} else {
		return "$marker_path/$bname.updated.ncbimap";
	}
}

=head2 get_count_from_reps

input marker name
Returns the number of representatives for a marker using the .rep file from the reference package
=cut

sub get_count_from_reps {
	my %args        = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_file = get_marker_rep_file( self => $self, marker => $marker );
	my $rep_num     = `grep -c '>' "$marker_file"`;
	chomp($rep_num);
	return $rep_num;
}

=head2 concatenate_alignments

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub concatenate_alignments {
	my %args          = @_;
	my $self          = $args{self};
	my $outputFasta   = $args{output_fasta};
	my $outputMrBayes = $args{output_bayes};
	my $gapmultiplier = $args{gap_multiplier};                           # 1 for protein, 3 for reverse-translated DNA
	my $aln_ref       = $args{alignments};
	my @alignments    = @$aln_ref;
	my $catobj        = 0;
	my $MRBAYES       = ps_open(">$outputMrBayes");
	my $partlist      = "partition genes = ".scalar(@alignments).": ";
	my $prevlen       = 0;
	foreach my $file (@alignments) {
		my $aln;
		my $marker = basename($file);
		$gapmultiplier = 1 if ( $marker =~ /16s/ || $marker =~ /18s/ );
		$marker =~ s/\..+//g;                                            # FIXME: this should really come from a list of markers
		unless ( -e $file ) {

			# this marker doesn't exist, need to create a dummy with the right number of gap columns
			$aln = make_dummy_file(
									self           => $self,
									marker         => $marker,
									gap_multiplier => $gapmultiplier
			);
		} else {
			my $in = Bio::AlignIO->new( -file => $file, '-format' => 'fasta' );
			unless ( $aln = $in->next_aln() ) {

				# empty marker alignment file, need to create a dummy with the right number of gap columns
				$aln = make_dummy_file(
										self           => $self,
										marker         => $marker,
										gap_multiplier => $gapmultiplier
				);
			}
		}
		my $csname = $file;
		$csname =~ s/\..+//g;
		$partlist .= "," if $catobj != 0;
		$partlist .= $csname;
		print $MRBAYES "charset $csname = ".( $prevlen + 1 )."-".( $prevlen + $aln->length() )."\n";
		$prevlen += $aln->length();
		my $prevseq = 0;
		my $newaln = $aln->select( 1, 1 );
		$newaln->verbose(-1);

		foreach my $curseq ( $aln->each_alphabetically() ) {
			if ( $prevseq == 0 ) {
				$prevseq = $curseq;
				next;
			}
			if ( $prevseq->id ne $curseq->id ) {
				$curseq->verbose(-1);
				$newaln->add_seq($curseq);
			}
			$prevseq = $curseq;
		}
		$aln = $newaln;
		if ( $catobj == 0 ) {
			$catobj = $aln;
			next;
		}

		# add any sequences missing from this sequence
		foreach my $catseq ( $catobj->each_alphabetically() ) {
			if ( length( $aln->each_seq_with_id( $catseq->id ) == 0 ) ) {

				# add this sequence as all gaps
				my $tmpseq = "-" x $aln->length();
				my $newseq = Bio::LocatableSeq->new(
													 -seq      => $tmpseq,
													 -alphabet => "protein",
													 -id       => $catseq->id,
													 start     => 0,
													 end       => 0
				);
				$newseq->verbose(-1);
				$aln->add_seq($newseq);
			}
		}

		# vice versa
		foreach my $alnseq ( $aln->each_alphabetically() ) {
			if ( length( $catobj->each_seq_with_id( $alnseq->id ) == 0 ) ) {

				# add this sequence as all gaps
				my $tmpseq = "-" x $catobj->length();
				my $newseq = Bio::LocatableSeq->new(
													 -seq      => $tmpseq,
													 -alphabet => "protein",
													 -id       => $alnseq->id,
													 start     => 0,
													 end       => 0
				);
				$newseq->verbose(-1);
				$catobj->add_seq($newseq);
			}
		}
		$catobj = cat( $catobj, $aln );
		$catobj->verbose(-1);
	}
	return unless $catobj;    # exit if there's nothing to write about
	print $MRBAYES "$partlist;\n";
	my $out = Bio::AlignIO->new( -file => ">$outputFasta", '-format' => 'fasta' );
	foreach my $dummyseq ( $catobj->each_seq_with_id("dummydummydummy") ) {
		$catobj->remove_seq($dummyseq);
	}
	$out->write_aln($catobj);

	# get rid of the trailing /1-blahblah that BioPerl insists on adding to sequence names
	# (naughty BioPerl! time out for you!!)
	my $FA       = ps_open($outputFasta);
	my @fa_lines = <$FA>;
	my $FA_REOUT = ps_open( ">".$outputFasta );
	foreach my $line (@fa_lines) {
		$line =~ s/\/\d-\d+//g;
		print $FA_REOUT $line;
	}
	close $FA_REOUT;
}

=head2 write_chunk_completion_to_run_info

Writes the completed chunk information to the run_info_file

=cut

sub write_step_completion_to_run_info {
	my %args    = @_;
	my $self    = $args{self} || miss("PS Object");
	my $chunk   = $args{chunk} || miss("Chunk");
	my $step    = $args{step} || miss("Step");
	my $RUNINFO = ps_open( ">>".get_run_info_file( self => $self ) );
	print $RUNINFO "Chunk $chunk $step completed\t"
	  .$self->{'run_info'}{$chunk}{$step}[0]."\t"
	  .$self->{'run_info'}{$chunk}{$step}[1]."\t"
	  .$self->{'run_info'}{$chunk}{$step}[2]."\n";
	close($RUNINFO);
}

=head2 start_step

Initializes the timer for a specific step for a specific chunk
Also adds the value into the run_info data structure

=cut

sub start_step {
	my %args              = @_;
	my $self              = $args{self} || miss("PS Object");
	my $step              = $args{step} || miss("Step to look for in run_info.txt\n");
	my $chunk             = $args{chunk} || miss("Chunk");
	my $start_search_time = start_timer( name => "start_$step"."_$chunk", silent => 1 );
	push( @{ $self->{'run_info'}{$chunk}{$step} }, $start_search_time );
}

=head2 end_step

Ends the timer for a specific step for a specific chunk
Also adds the values into the run_info data structure

=cut

sub end_step {
	my %args            = @_;
	my $self            = $args{self} || miss("PS Object");
	my $step            = $args{step} || miss("Step to look for in run_info.txt\n");
	my $chunk           = $args{chunk} || miss("Chunk");
	my $end_search_time = start_timer( name => "end_$step"."_$chunk", silent => 1 );
	push( @{ $self->{'run_info'}{$chunk}{$step} }, $end_search_time, end_timer( name => "start_$step"."_$chunk", silent => 1 ) );
}

=head2 load_run_info

Reads the run_info.txt file from a run and loads it into a data structure

=cut

sub load_run_info {
	my %args        = @_;
	my $self        = $args{self};
	my %return_hash = ();
	return \%return_hash unless -e get_run_info_file( self => $self );
	my $RUN_INFO_FH = ps_open( get_run_info_file( self => $self ) );
	while (<$RUN_INFO_FH>) {
		chomp($_);

		#debug $_."\n";
		next unless $_ =~ m/^Chunk.*completed/;
		my @line = split( /\t/, $_ );
		## data structure run_info{chunk}{}
		$line[0] =~ m/Chunk\s+(\d+)\s+(\S+)\s+/;
		my $chunk    = $1;
		my $mode     = $2;
		my $start    = $line[1];
		my $end      = $line[2];
		my $duration = $line[3];
		push( @{ $return_hash{$chunk}{$mode} }, $start, $end, $duration );
	}
	return \%return_hash;
}

=head2 has_step_completed

Checks to see if a chunk has completed the step specified
returns 1 if it has and 0 if it hasn't
Always return 0 if $Phylosift::Settings::force is set to 1
=cut

sub has_step_completed {
	my %args  = @_;
	my $self  = $args{self} || miss("PS object");
	my $step  = $args{step} || miss("Step to look for in run_info.txt\n");
	my $chunk = $args{chunk} || miss("Chunk");
	my $force = $args{force};
	if ( defined( $self->{'run_info'}{$chunk}{$step} ) && scalar( @{ $self->{'run_info'}{$chunk}{$step} } ) > 1 ) {
		cleanup_chunk( self => $self, step => $step, chunk => $chunk ) if $force;
		return 1 unless $force;
	}
	return 0;
}

=head2 cleanup_chunk

Removes chunk data from a specific step

=cut

sub cleanup_chunk {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $step       = $args{step} || miss("Step to clean up");
	my $chunk      = $args{chunk} || miss("Chunk targetted");
	my $marker_ref = $args{markers};
	my @files_list = ();
	my $file;
	if ( $step eq 'Search' ) {
		push( @files_list, get_search_output_all_candidate( self => $self, chunk => $chunk ) );
	} elsif ( $step eq 'Align' ) {
		push( @files_list, get_aligner_output_all_unmasked( self => $self, chunk => $chunk ) );
		push( @files_list, get_aligner_output_all_newCandidate( self => $self, chunk => $chunk ) );
		push( @files_list, get_aligner_output_all_fasta( self => $self, chunk => $chunk ) );
	} elsif ( $step eq 'Place' ) {
		push( @files_list, get_place_output_all_jplace( self => $self, chunk => $chunk ) );
	} elsif ( $step eq 'Summarize' ) {
		push( @files_list, get_taxa_90pct_HPD( self => $self ) );
		push( @files_list, get_taxasummary( self => $self ) );
		push( @files_list, get_summarize_output_sequence_taxa( self => $self, chunk => $chunk ) );
		push( @files_list, get_summarize_output_sequence_taxa_summary( self => $self, chunk => $chunk ) );
	}
	return if @files_list == 0;    # no need to do anything else if the list is empty
	my $cmd = "rm ".join( ' ', @files_list );
	`$cmd`;
}

=head2 get_search_output_all_candidate
return an array of all candidate files in the blastDir for a specific chunk
=cut

sub get_search_output_all_candidate {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $chunk      = $args{chunk} || miss("Chunk");
	my @files_list = ();
	push( @files_list, glob( $self->{'blastDir'}."/*.candidate.*.$chunk.*" ) );
	return @files_list;
}

=head2 get_aligner_output_all_fasta
return an array of all fasta files in the alignDir for a specific chunk
=cut

sub get_aligner_output_all_fasta {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $chunk      = $args{chunk} || miss("Chunk");
	my @files_list = ();
	push( @files_list, glob( $self->{'alignDir'}."/*.$chunk.fasta" ) );
	return @files_list;
}

=head2 get_aligner_output_all_newCandidate
return an array of all newCandidate files in the alignDir
=cut

sub get_aligner_output_all_newCandidate {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $chunk      = $args{chunk} || miss("Chunk");
	my @files_list = ();
	push( @files_list, glob( $self->{'alignDir'}."/*.newCandidate.aa.$chunk" ) );
	return @files_list;
}

=head2 get_aligner_output_all_unmasked
return an array of all unmasked files in the alignDir
=cut

sub get_aligner_output_all_unmasked {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $chunk      = $args{chunk} || miss("Chunk");
	my @files_list = ();
	push( @files_list, glob( $self->{'alignDir'}."/*.$chunk.unmasked" ) );
	return @files_list;
}

=head2 get_place_output_all_jplace
return an array of all jplace files in the treeDir for a specific chunk
=cut

sub get_place_output_all_jplace {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $chunk      = $args{chunk} || miss("Chunk");
	my @files_list = ();
	push( @files_list, glob( $self->{'treeDir'}."/*.$chunk.jplace" ) );
	return @files_list;
}

=head2 get_summarize_output_sequence_taxa
return an array of all sequence_taxa files in the workingDir for a specific chunk
if chunk isn't defined, return all sequence_taxa files
=cut

sub get_summarize_output_sequence_taxa {
	my %args  = @_;
	my $self  = $args{self} || miss("PS object");
	my $chunk = $args{chunk};
	$chunk = '*' unless ( defined $chunk );
	my @files_list = ();
	push( @files_list, glob( $Phylosift::Settings::file_dir."/sequence_taxa.$chunk.txt" ) );
	return @files_list;
}

=head2 get_summarize_output_sequence_taxa_summary
return an array of all sequence_taxa_summary files in the workingDir for a specific chunk
if chunk isn't defined, return all sequence_taxa_summary files
=cut

sub get_summarize_output_sequence_taxa_summary {
	my %args  = @_;
	my $self  = $args{self} || miss("PS object");
	my $chunk = $args{chunk};
	$chunk = '*' unless ( defined $chunk );
	my @files_list = ();
	push( @files_list, glob( $Phylosift::Settings::file_dir."/sequence_taxa_summary.$chunk.txt" ) );
	return @files_list;
}

=head2 get_taxasummary
return the path to taxasummary for the PS object
=cut

sub get_taxasummary {
	my %args = @_;
	my $self = $args{self} || miss("PS object");
	return $Phylosift::Settings::file_dir."/taxasummary.txt";
}

=head2 get_taxa_90pct_HPD
return the path to taxa_90pct_HPD for the PS object
=cut

sub get_taxa_90pct_HPD {
	my %args = @_;
	my $self = $args{self} || miss("PS object");
	return $Phylosift::Settings::file_dir."/taxa_90pct_HPD.txt";
}

=head2 concatenate_summary_files

Concatenates a set of files not includin the header line for the first file only

=cut

sub concatenate_summary_files {
	my %args            = @_;
	my $self            = $args{self} || miss("PS Object");
	my $output_file     = $args{output_file} || miss("Output file");
	my $files_array_ref = $args{files};
	my $OUTPUT_FH       = ps_open(">$output_file");
	my @files           = @{$files_array_ref};
	my $header          = "";
	my $printed_header  = 0;
	foreach my $file (@files) {
		my $FH = ps_open($file);
		$header = <$FH>;
		print $OUTPUT_FH $header unless $printed_header;
		$printed_header = 1;
		while (<$FH>) {    #print the rest of the file
			print $OUTPUT_FH $_;
		}
		close($FH);
	}
	close($OUTPUT_FH);
}

sub start_timer {
	my %args      = @_;
	my $timername = $args{name};
	my $silent    = $args{silent} || 0;
	my $t         = time;
	my @timerval  = localtime($t);
	$timers{$timername} = $t;
	my $readable_time =
	  sprintf( "%4d-%02d-%02d %02d:%02d:%02d", $timerval[5] + 1900, $timerval[4] + 1, $timerval[3], $timerval[2], $timerval[1], $timerval[0] );
	debug "Before $timername $readable_time\n" unless $silent;
	return $readable_time;
}

sub end_timer {
	my %args      = @_;
	my $timername = $args{name};
	my $silent    = $args{silent} || 0;
	my $t         = time;
	my @timerval  = localtime($t);
	debug sprintf( "After $timername %4d-%02d-%02d %02d:%02d:%02d\n",
				   $timerval[5] + 1900,
				   $timerval[4] + 1,
				   $timerval[3], $timerval[2], $timerval[1], $timerval[0] )
	  unless $silent;
	return $t - $timers{$timername};
}

=head2 marker_oldstyle

Checks whether a marker is in the old style format (PMPROK*) or the new format (marker package directory)

=cut

sub marker_oldstyle {
	my %args   = @_;
	my $marker = $args{markers};
	return 1 if ( $marker =~ /PMPROK/ );
	return 0;
}

=head2 escape_char

Inserts \ before problematic characters 
Characters currently flagged as problematic : < > (space)
Add as needed.

=cut

sub escape_char {
	my %args   = @_;
	my $string = $args{string};
	$string =~ s/\s/\\ /g;
	return $string;
}

=head2 open_SeqIO_object

Opens a sequence file and returns a SeqIO object.  Allows for gzip and bzip compression
returns a SeqIO object

=cut

sub open_SeqIO_object {
	my %args   = @_;
	my $format = $args{format} || "FASTA";      #default
	my $file   = $args{file} || miss("file");
	my $io_object;
	if ( exists $args{format} ) {
		$format = $args{format};
	}
	if ( $args{file} =~ /\.gz$/ ) {
		$io_object = Bio::SeqIO->new( -file   => "gzip -cd \"$args{file}\" |",
									  -format => $format );
	} elsif ( $args{file} =~ /\.bz2$/ ) {
		$io_object = Bio::SeqIO->new( -file   => "bzcat \"$args{file}\" |",
									  -format => $format );
	} else {
		$io_object = Bio::SeqIO->new( -file => $args{file}, -format => $format );
	}
	return $io_object;
}

=head2 open_sequence_file

Opens a sequence file, either directly or by decompressing it with gzip or bzip2
Returns an open filehandle

=cut

sub open_sequence_file {
	my %args = @_;
	my $file = $args{file} || miss("file");
	my $F1IN;
	if ( $file =~ /\.gz$/ ) {
		$F1IN = ps_open("gzip -cd \"$file\" |");
	} elsif ( $file =~ /\.bz2$/ ) {
		$F1IN = ps_open("bzcat \"$file\" |");
	} elsif ( $file eq "STDIN" ) {
		$F1IN = *STDIN;
	} else {
		$F1IN = ps_open($file);
	}
	return $F1IN;
}

sub get_db {
	my %args    = @_;
	my $path    = $args{path} || $Phylosift::Settings::marker_dir;
	my $self    = $args{self};                                       # optional argument;
	my $db_name = $args{db_name};
	if ( defined($self) && !defined( $args{path} ) ) {
		return $Phylosift::Settings::markers_extended_dir."/$db_name"
		  if defined($Phylosift::Settings::extended) && $Phylosift::Settings::extended;
		return $Phylosift::Settings::marker_dir."/$db_name";
	}
	return "$path/$db_name";
}

=head2 read_interleaved_fastq

Reads 1 fastq file with interleaved paired ends reads
Writes to the provided pipe in fasta format.

=cut

sub read_interleaved_fastq {
	my %args       = @_;
	my $pipe       = $args{PIPEOUT};
	my $suffix_1   = $args{suffix_1};
	my $suffix_2   = $args{suffix_2};
	my $PIPE_IN    = open_sequence_file( file => $args{file} );
	my $read_id_1  = <$PIPE_IN>;
	my $read_seq_1 = <$PIPE_IN>;
	my $read_id_2  = <$PIPE_IN>;
	my $read_seq_2 = <$PIPE_IN>;
	$read_id_1 =~ m/^\@(\S+)$suffix_1$/;
	my $core_1 = $1;
	$read_id_2 =~ m/^\@(\S+)$suffix_2$/;
	my $core_2 = $1;

	while (<$PIPE_IN>) {
	}
}

=head2 get_blastp_db

Returns the name and path of the blast DB

=cut

sub get_blastp_db {
	my %args = @_;
	$args{db_name} = "blastrep";
	return get_db(%args);
}

=head2 get_rna_db

returns the name and path of the RNA DB

=cut

sub get_rna_db {
	my %args = @_;
	$args{db_name} = "rnadb";
	return get_db(%args);
}

=head2 get_lastal_db

returns the name and path of the lastal DB

=cut

sub get_lastal_db {
	my %args = @_;
	$args{db_name} = "replast";
	return get_db(%args);
}

sub get_candidate_file {
	my %args   = @_;
	my $self   = $args{self};
	my $marker = $args{marker} || miss("Marker name");
	my $type   = $args{type};
	my $dna    = $args{dna};
	my $new    = $args{new};
	my $chunk  = $args{chunk};
	my $ffn    = ".aa";
	$ffn = ".ffn" if ( defined( $args{dna} ) && $dna );
	my $candidate = ".candidate";
	$candidate = ".newCandidate" if defined( $args{new} ) && $new;
	my $dir = $self->{"blastDir"};
	$dir = $self->{"alignDir"} if defined( $args{new} ) && $new;
	my $chunky = "";
	$chunky = ".$chunk" if defined($chunk);
	$marker =~ s/.+\///g;    # strip off any prepended directories
	return "$dir/$marker$type$candidate$ffn$chunky";
}

=head2 index_marker_db

Indexes the marker database for searches with Lastal
Input: marker list and self

=cut

sub index_marker_db {
	my %args    = @_;
	my $self    = $args{self};
	my $path    = $args{path};
	my @markers = @{ $args{markers} };
	debug "Indexing ".scalar(@markers)." in $path\n";
	debug "Inside index_marker_db\n";

	# use alignments to make an unaligned fasta database containing everything
	# strip gaps from the alignments
	my $rna_db        = get_rna_db( path => $path );
	my $rna_db_fasta  = "$rna_db.fasta";
	my $PDBOUT        = ps_open( ">".get_blastp_db( path => $path ) );
	my $RNADBOUT      = ps_open( ">".$rna_db_fasta );
	my $MARKERLISTOUT = ps_open(">$path/marker_list.txt");
	foreach my $marker (@markers) {
		my $marker_rep = get_marker_rep_file(
											  self    => $args{self},
											  marker  => $marker,
											  updated => 1
		);
		$marker_rep = get_marker_rep_file(
										   self    => $args{self},
										   marker  => $marker,
										   updated => 0
		) unless -e $marker_rep;

		#debug "marker $marker is protein\n"
		#if is_protein_marker( marker => $marker );
		#debug "marker rep file $marker_rep\n";
		my $DBOUT = $RNADBOUT;
		$DBOUT = $PDBOUT if is_protein_marker( marker => $marker );
		next if ( !is_protein_marker( marker => $marker ) && $marker =~ /\.short/ );
		unless ( -f $marker_rep ) {
			warn "$marker_rep not found.\n" unless $marker_rep =~ m/1[68]s_reps(_\w\w\w)?.short./;
			warn "Warning: marker $marker appears to be missing data\n" if ( $marker !~ /\.short/ );
			next;
		}
		print $MARKERLISTOUT "$marker\n";
		my $INALN = ps_open($marker_rep);
		while ( my $line = <$INALN> ) {
			if ( $line =~ /^>(.+)/ ) {
				print $DBOUT "\n>$marker"."__$1\n";
			} else {
				$line =~ s/[-\.\n\r]//g;
				$line =~ tr/a-z/A-Z/;
				print $DBOUT $line;
			}
		}
	}
	print $PDBOUT "\n";
	print $RNADBOUT "\n";
	close $MARKERLISTOUT;
	close $PDBOUT;    # be sure to flush I/O
	close $RNADBOUT;
	my $blastp_db = get_blastp_db( path => $path );
	`mv "$blastp_db" "$path/rep.dbfasta"`;

	# make a last database
	`cd "$path" ; $Phylosift::Settings::lastdb -s 900M -p replast rep.dbfasta`;
	unlink("$path/rep.dbfasta");    # don't need this anymore!

	# make a rna database
	if ( -e $rna_db_fasta && -s $rna_db_fasta > 100 ) {
		`cd "$path" ; $Phylosift::Settings::lastdb "$rna_db" "$rna_db_fasta"`;
	}

	# now create the .hmm files if they aren't already present
	# this is the case in the extended marker set, since the hmms are too big for transit
	foreach my $marker (@markers) {
		my $hmm_file = get_marker_hmm_file( self => $args{self}, marker => $marker );
		next if -e $hmm_file;
		build_hmm( marker => $marker ) unless $marker =~ /1[68]s_reps(_\w\w\w)?.short/;
	}
	return 1;
}

sub build_hmm {
	my %args     = @_;
	my $marker   = $args{marker} || miss("marker");
	my $hmm_file = get_marker_hmm_file( self => $args{self}, marker => $marker );
	my $stk_file = get_marker_stockholm_file( self => $args{self}, marker => $marker );
	`$Phylosift::Settings::hmmbuild "$hmm_file" "$stk_file"`;
}

=head2 gather_markers
=over

=item *

    ARGS : $markerFile - file to read the marker names from that will be used in the pipeline.
    Reads markers from a file if specified by the user OR gathers all the markers from the markers directory
    The Marker names are stored in an Array
    If the filename is empty, use all the default markers.

=back

=cut

sub gather_markers {
	my %args         = @_;
	my $marker_file  = $args{marker_file};
	my $missing_hmm  = $args{allow_missing_hmm} || 0;
	my $force_gather = $args{force_gather} || 0;
	my $path         = $args{path} || $Phylosift::Settings::marker_dir;
	my @marks        = ();

	# try to use a marker list file if it exists
	if ( !defined($marker_file) && !$force_gather ) {
		$marker_file = "$path/marker_list.txt" if -e "$path/marker_list.txt";
	}

	#create a file with a list of markers called markers_list.txt
	if ( defined($marker_file) && -f $marker_file && $marker_file ne "" ) {
		debug "Using a marker list file $marker_file\n";

		#gather a custom list of makers, list convention is 1 marker per line
		my $MARKERS_IN = ps_open($marker_file);
		while (<$MARKERS_IN>) {
			chomp($_);
			push( @marks, $_ );
			$marker_lookup{ get_marker_basename( marker => $_ ) } = $_;
		}
		close($MARKERS_IN);
		return @marks;
	} else {

		# gather all markers
		# this is for the original marker set
		my @files = <$path/*.faa>;
		foreach my $file (@files) {
			$file =~ m/\/(\w+).faa/;
			$marker_lookup{ get_marker_basename( marker => $1 ) } = $1;
			push( @marks, $1 );
		}

		# now gather directory packaged markers (new style)
		# use maxdepth 2 for two-directory-level markers in the extended marker set
		my $MLIST = ps_open("find $path -maxdepth 2 -mindepth 1 -type d |");
		while ( my $line = <$MLIST> ) {
			chomp $line;

			#debug "$line\n";
			next if $line =~ /PMPROK/;
			next if $line =~ /concat/;
			next if $line =~ /representatives/;
			next if $line =~ /.updated$/;                # just include the base version name
			next if $line =~ /codon.updated.sub\d+$/;    # just include the base version name
			$line = substr( $line, length($path) + 1 );

			# all markers need to have an hmm or a cm else they are not usable
			my $baseline = $line;
			$baseline =~ s/.+\///g;
			if ( !$missing_hmm ) {
				next unless ( -e "$path/$line/$line.cm" || -e "$path/$line/$baseline.hmm" );
			}
			$marker_lookup{$baseline} = $line;
			push( @marks, $line );
		}
		debug( "Found ".scalar(@marks)." markers\n" );
	}

	# always sort these, since the concats need to be used in the same order every time
	return sort @marks;
}

=head2 get_available_memory

Reads total installed memory on the system

=cut

sub get_available_memory {
	my $mem = 4000000;    #default 4Gigs
	if ( $^O =~ m/darwin/ ) {
		my $inf = `/usr/sbin/system_profiler SPHardwareDataType | grep Memory`;
		if ( $inf =~ /      Memory: (\d+) GB/ ) {
			$mem = $1 * 1048576;
		}
	} else {
		my $inf = `cat /proc/meminfo | grep MemTotal`;
		if ( $inf =~ /MemTotal:\s+(\d+) kB/ ) {
			$mem = $1;
		}
	}
	return $mem;
}

=head2 get_sequence_input_type

Checks whether input is FastA, FastQ, which quality type (33 or 64), and DNA or AA
Returns a hash reference with the values 'seqtype', 'format', and 'qtype' populated.

=cut

sub get_sequence_input_type {
	my $FILE = shift;
	my %type;
	my $counter    = 0;
	my $maxfound   = 0;
	my $dnacount   = 0;
	my $line_count = 0;
	$type{seqtype} = "dna";
	$type{format}  = "unknown";
	$type{qtype}   = "none";
	my $allcount = 0;
	my $sequence = 1;
	my $minq     = 255;    # minimum fastq quality score (for detecting phred33/phred64)
	my @lines;

	while ( my $line = <$FILE> ) {
		push( @lines, $line );
		if ( $line =~ /^>/ ) {
			$maxfound = $counter > $maxfound ? $counter : $maxfound;
			$counter = 0;
			$type{format} = "fasta" if $type{format} eq "unknown";
		} elsif ( $line =~ /^@/ || $line =~ /^\+/ ) {
			$counter  = 0;
			$sequence = 1;
			$type{format} = "fastq" if $type{format} eq "unknown";
		} elsif ( $line =~ /^\+/ ) {
			$sequence = 0;
			$type{format} = "fastq" if $type{format} eq "unknown";
		} elsif ( $type{format} eq "fastq" && !$sequence ) {

			# check whether qualities are phred33 or phred64
			for my $q ( split( //, $line ) ) {
				$minq = ord($q) if ord($q) < $minq;
			}
		} elsif ($sequence) {

			#$sequence =~ s/[-\.]//g;    #removing gaps from the sequences
			$line =~ s/[-\.]//g;
			$counter  += length($line) - 1;
			$dnacount += $line =~ tr/[ACGTUNacgtun]//;
			$allcount += length($line) - 1;
		}
		$line_count++;
		last if ( $line_count > 1000 );
		last if ( $counter > 100000 );
	}

	$maxfound = $counter > $maxfound ? $counter : $maxfound;
	$type{seqtype} = "protein" if ( $dnacount < $allcount * 0.75 );
	$type{seqtype} = "dna"
	  if ( $type{format} eq "fastq" );    # nobody using protein fastq (yet)
	$type{qtype} = "phred64" if $minq < 255;
	$type{qtype} = "phred33" if $minq < 64;
	$type{paired} = 0;                    # TODO: detect interleaved read pairing
	$type{buffer} = \@lines;              # a copy of all lines read in case we're streaming
	return \%type;
}

=head2 get_sequence_input_type_quickndirty

Checks whether input is either short sequence reads, e.g. < 500nt or assembled fragments
without reading the whole file.

=cut

sub get_sequence_input_type_quickndirty {
	my $file         = shift;
	my $maxshortread = 500;
	my $FILE         = ps_open($file);
	my $filesize     = -s "$file";
	my $counter      = 0;
	my $maxfound     = 0;
	my $allcount     = 0;
	my $dnacount     = 0;
	my $seqtype      = "dna";
	my $length       = "long";
	my $format       = "unknown";
	for ( my $i = 0; $i < 200; $i++ ) {
		my $seekpos = int( rand( $filesize - 100 ) );
		$seekpos = 0
		  if ( $i == 0 );    # always start with the first line in case the sequence is on a single line!
		seek( $FILE, $seekpos, 0 );
		$counter = 0;
		my $line = <$FILE>;    # burn a line to be sure we get to sequence
		while ( $line = <$FILE> ) {
			if ( $line =~ /^>/ ) {
				$format = "fasta";
				last if $i > 0;
				$i++;
			} elsif ( $line =~ /^@/ || $line =~ /^\+/ ) {
				$format = "fastq";
				last if $i > 0;
				$i++;
			} else {
				$counter += length($line);
				$dnacount += $line =~ tr/[ACGTNacgtn]//;
			}
		}
		$maxfound = $counter if $maxfound < $counter;
		$allcount += $counter;
		last if ( $counter > 500 );    # found a long read
	}
	$seqtype = "protein" if ( $dnacount < $allcount * 0.75 );
	$seqtype = "dna"
	  if ( $format eq "fastq" );       # nobody using protein fastq (yet)
	my $aamult = $seqtype eq "protein" ? 3 : 1;
	$length = "short" if $maxfound < ( 500 / $aamult );
	return ( $seqtype, $length, $format );
}

=head2 get_date_YYYYMMDD

Gets the current date and formats it by YYYYMMDD

=cut

sub get_date_YYYYMMDD {
	my @timerval = localtime();
	my $datestr  = ( 1900 + $timerval[5] );
	$datestr .= 0 if $timerval[4] < 9;
	$datestr .= ( $timerval[4] + 1 );
	$datestr .= 0 if $timerval[3] <= 9;
	$datestr .= $timerval[3];
	return $datestr;
}

=head2 print_citations

prints out suggested citations for the analysis

=cut

sub print_citations {
	print "PhyloSift -- Phylogenetic analysis of genomes and metagenomes\n";
	print "(c) 2011-2014 Aaron Darling and Guillaume Jospin\n";
	print "\nCITATION:\n";
	print "		Darling AE, Jospin G, Lowe E, Matsen FA IV, Bik HM, Eisen JA (2014)\n";
	print "		PhyloSift: phylogenetic analysis of genomes and metagenomes. PeerJ 2:e243 http://dx.doi.org/10.7717/peerj.243\n";
	print "\n\nPhyloSift incorporates several other software packages, please consider also citing the following papers:\n";
	print qq{

		pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree.
		Frederick A Matsen, Robin B Kodner, and E Virginia Armbrust
		BMC Bioinformatics 2010, 11:538
		
		Adaptive seeds tame genomic sequence comparison.
		SM Kielbasa, R Wan, K Sato, P Horton, MC Frith
		Genome Research 2011.
		
		Infernal 1.0: Inference of RNA alignments
		E. P. Nawrocki, D. L. Kolbe, and S. R. Eddy
		Bioinformatics 25:1335-1337 (2009)
 		
 		Bowtie: Ultrafast and memory-efficient alignment of short DNA sequences to the human genome.
 		Langmead B, Trapnell C, Pop M, Salzberg SL. Genome Biol 10:R25.
 		
 		HMMER 3.0 (March 2010); http://hmmer.org/
 		Copyright (C) 2010 Howard Hughes Medical Institute.
 		Freely distributed under the GNU General Public License (GPLv3).
 		
 		Phylogenetic Diversity within Seconds.
		Bui Quang Minh, Steffen Klaere and Arndt von Haeseler
		Syst Biol (2006) 55 (5): 769-773.
};
}

=head2 alignment_to_fasta

intput: alignment file , output path
Removes all gaps from an alignment file and writes the sequences in fasta format in the target directory

=cut

sub alignment_to_fasta {
	my %args       = @_;
	my $aln_file   = $args{aln};
	my $target_dir = $args{target_dir};
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $in = open_SeqIO_object( file => $aln_file );
	my $FILEOUT = ps_open( ">".$target_dir."/".$core.".fasta" );
	while ( my $seq_object = $in->next_seq() ) {
		my $seq = $seq_object->seq;
		my $id  = $seq_object->id;
		$seq =~ s/-\.//g;    # shouldnt be any gaps
		print $FILEOUT ">".$id."\n".$seq."\n";
	}
	close($FILEOUT);
	return "$target_dir/$core.fasta";
}

=head2 generate_hmm

input: alignment_file, target_directory
generates a HMM profile from an alignement in FASTA format (arg) using hmmbuild.  The hmm is placed in the target_directory

=cut

sub generate_hmm {
	my %args       = @_;
	my $file_name  = $args{file};
	my $target_dir = $args{target_dir};
	debug("HMM on : $file_name\n");
	my ( $core_name, $path, $ext ) = fileparse( $file_name, qr/\.[^.]*$/ );
	if ( !-e "$target_dir/$core_name.hmm" ) {
		`hmmbuild --informat afa "$target_dir/$core_name.hmm" "$file_name"`;
	}
	return "$target_dir/$core_name.hmm";
}

=head2 hmmalign_to_model

input : hmm_profile,sequence_file,target_dir
Aligns sequences to an HMM model and outputs an alignment
=cut

sub hmmalign_to_model {
	my %args          = @_;
	my $hmm_profile   = $args{hmm};
	my $sequence_file = $args{sequence_file};
	my $target_dir    = $args{target_dir};
	my $ref_ali       = $args{reference_aln};
	my ( $core_name, $path, $ext ) = fileparse( $sequence_file, qr/\.[^.]*$/ );
	if ( !-e "$target_dir/$core_name.aln" ) {
		`hmmalign --mapali "$ref_ali" --trim --outformat afa -o "$target_dir/$core_name.aln" "$hmm_profile" "$sequence_file"`;

		#	    `hmmalign --outformat afa -o $target_dir/$core_name.aln $hmm_profile $sequence_file`;
	}
	return "$target_dir/$core_name.aln";
}

=head2 unalign_sequences

input : aln_file
Masks the unaligned columns out of an alignemnt file. Removes ( and ) from the sequence names 
Also removes duplicate IDs
=cut

sub unalign_sequences {
	my %args        = @_;
	my $aln_file    = $args{aln};
	my $output_path = $args{output_path};
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $in        = open_SeqIO_object( file => $aln_file );
	my $seq_count = 0;
	my $FILEOUT   = ps_open(">$output_path");
	while ( my $seq_object = $in->next_seq() ) {
		my $seq = $seq_object->seq;
		my $id  = $seq_object->id;
		$seq =~ s/[-\.]//g;    # shouldnt be any gaps
		print $FILEOUT ">".$id."\n".$seq."\n";
		$seq_count++;
	}
	close($FILEOUT);
	return $seq_count;
}

=head2 get_run_info_file

returns the path to the run_info_file

=cut

sub get_run_info_file {
	my %args = @_;
	my $self = $args{self} || miss("PS object");
	return $Phylosift::Settings::file_dir."/run_info.txt";
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

    perldoc Phylosift::Utilities


You can also look for information at:

=over 4

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

1;    # End of Phylosift::Utilities
