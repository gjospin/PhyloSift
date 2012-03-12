package Phylosift::Utilities;

#use 5.006;
use strict;
use warnings;
use FindBin qw($Bin);
BEGIN { unshift( @INC, "$FindBin::Bin/../legacy/" ) if $] < 5.01; }
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Align::Utilities qw(:all);
use Bio::TreeIO;
use Bio::Tree::Tree;
use POSIX ();
use LWP::Simple;
use Carp;
use Cwd;
require File::Fetch;

if ( $^O =~ /arwin/ ) {
	use lib "$FindBin::Bin/../osx/darwin-thread-multi-2level/";
}
use Exporter;
use vars qw[ @EXPORT @EXPORT_OK %EXPORT_TAGS @ISA ];
@ISA       = 'Exporter';
@EXPORT    = qw[start_timer end_timer debug];
@EXPORT_OK = qw[];
%EXPORT_TAGS = (
				 STD => \@EXPORT,
				 all => [ @EXPORT, @EXPORT_OK ],
);
our $debuglevel = 0;

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

=head2 programChecks

checks the program requirements for PhyloSift
writes to STDERR if a program is missing or the wrong version is installed
returns 1 or 0 depending on success of failure.

=cut

sub get_program_path {
	my %args =@_;
	my $progname  = $args{prog_name};
	my $progpath  = $args{prog_path};
	my $progcheck = "";

	#	print STDERR "BIN : $Bin\n";
	#	exit;
	if ( defined($progpath) && $progpath ne "" && -x $progpath . "/" . $progname ) {
		$progcheck = $progpath . "/" . $progname;
	} else {
		$progcheck = `which $progname 2> /dev/null`;
		chomp $progcheck;
	}

	# last ditch attempt, check the directories from where the script is running
	$progcheck = $Bin . "/" . $progname     unless ( $progcheck =~ /$progname/ || !( -x $Bin . "/" . $progname ) );
	$progcheck = $Bin . "/bin/" . $progname unless ( $progcheck =~ /$progname/ || !( -x $Bin . "/bin/" . $progname ) );

	# check the OS and use Mac binaries if needed
	if ( $^O =~ /arwin/ ) {
		$progcheck = $Bin . "/osx/" . $progname unless ( $progcheck =~ /$progname/ && !( -x $Bin . "/" . $progname ) );
		$progcheck = $Bin . "/osx/" . $progname if ( $progcheck =~ /$Bin\/bin/ );    # don't use the linux binary!
	}
	return $progcheck;
}

# external programs used by Phylosift
our $pplacer      = "";
our $guppy        = "";
our $rppr         = "";
our $taxit        = "";
our $hmmalign     = "";
our $hmmsearch    = "";
our $hmmbuild     = "";
our $blastall     = "";
our $formatdb     = "";
our $rapSearch    = "";
our $preRapSearch = "";
our $raxml        = "";
our $readconciler = "";
our $bowtie2align = "";
our $bowtie2build = "";
our $cmalign      = "";
our $pda          = "";
our $fasttree     = "";
our $lastdb       = "";
our $lastal       = "";

sub programChecks {
	eval 'require Bio::Seq;';
	if ($@) {
		carp "Bioperl was NOT found\n";
		return 1;
	}
	$pplacer = get_program_path( prog_name=>"pplacer", prog_path=>$Phylosift::Settings::pplacer_path );
	if ( $pplacer eq "" ) {

		#program not found return;
		carp("pplacer not found");
		return 1;
	} else {
		`$pplacer --version` =~ m/v1.1.alpha(\d+)/;
		if ( $1 < 10 ) {

			# pplacer was found but the version doesn't match the one tested with Phylosift
			carp("Warning : a different version of pplacer was found. PhyloSift was tested with pplacer v1.1.alpha10\n");
		}
	}
	$guppy    = get_program_path( prog_name=>"guppy",   prog_path=> $Phylosift::Settings::ps_path );
	$taxit    = get_program_path( prog_name=>"taxit",    prog_path=>$Phylosift::Settings::ps_path );
	$rppr     = get_program_path( prog_name=>"rppr",    prog_path=> $Phylosift::Settings::ps_path );
	$cmalign  = get_program_path( prog_name=>"cmalign", prog_path=> $Phylosift::Settings::ps_path );
	$hmmalign = get_program_path( prog_name=>"hmmalign",prog_path=> $Phylosift::Settings::hmmer3_path );
	if ( $hmmalign eq "" ) {

		#program not found return;
		carp("HMMER3 not found");
		return 1;
	} elsif ( `$hmmalign -h` !~ m/HMMER 3.0rc1/ ) {

		# pplacer was found but the version doens't match the one tested with Phylosift
		carp "Warning : a different version of HMMER was found. PhyloSift was tested with HMMER 3.0rc1\n";
	}
	$hmmsearch = get_program_path( prog_name=>"hmmsearch",prog_path=> $Phylosift::Settings::hmmer3_path );
	$hmmbuild  = get_program_path( prog_name=>"hmmbuild", prog_path=> $Phylosift::Settings::hmmer3_path );
	$rapSearch = get_program_path( prog_name=>"rapsearch",prog_path=> $Phylosift::Settings::ps_path );
	if ( $rapSearch eq "" ) {
		carp("rapsearch was not found\n");
		return 1;
	}
	$preRapSearch = get_program_path( prog_name=>"prerapsearch", prog_path=>$Phylosift::Settings::ps_path );
	$blastall     = get_program_path( prog_name=>"blastall",   prog_path=>  $Phylosift::Settings::ps_path );
	if ( $blastall eq "" ) {
		carp("blastall was not found\n");
		return 1;
	}
	$formatdb = get_program_path( prog_name=>"formatdb",prog_path=> $Phylosift::Settings::ps_path );
	$raxml    = get_program_path( prog_name=>"raxmlHPC", prog_path=>$Phylosift::Settings::ps_path );
	if ( $raxml eq "" ) {
		carp("raxmlHPC was not found\n");
		return 1;
	}
	$readconciler = get_program_path( prog_name=>"readconciler", prog_path=> $Phylosift::Settings::ps_path );
	$pda          = get_program_path( prog_name=>"pda",          prog_path=> $Phylosift::Settings::ps_path );
	$fasttree     = get_program_path( prog_name=>"FastTree",     prog_path=> $Phylosift::Settings::ps_path );
	$lastdb       = get_program_path( prog_name=>"lastdb",       prog_path=> $Phylosift::Settings::ps_path );
	$lastal       = get_program_path( prog_name=>"lastal",       prog_path=> $Phylosift::Settings::ps_path );
	$bowtie2align = get_program_path( prog_name=>"bowtie2-align",prog_path=> $Phylosift::Settings::bowtie2_path );
	if ( $bowtie2align eq "" ) {

		#program not found return;
		carp("bowtie2 not found");
		return 1;
	}
	$bowtie2build = get_program_path(prog_name=> "bowtie2-build",prog_path=> $Phylosift::Settings::bowtie2_path );
	return 0;
}

=head2 dataChecks

Check for requisite PhyloSift marker datasets

=cut

our $marker_dir           = "";
our $markers_extended_dir = "";
our $ncbi_dir             = "";

sub get_data_path {
	my %args = @_;
	my $dataname  = $args{data_name};
	my $datapath  = $args{data_path};
	my $datacheck = "";
	if ( defined($datapath) && $datapath ne "" ) {
		$datacheck = $datapath . "/" . $dataname;
	} else {
		my $scriptpath = dirname($0);
		$scriptpath =~ s/bin\/?$//g;
		$datacheck = $scriptpath . "/share/phylosift/" . $dataname;
		return $datacheck if ( -x $datacheck );

		# if the system data dir doesn't exist, default to the user's home
		$datacheck = $ENV{"HOME"} . "/share/phylosift/" . $dataname;
	}
	return $datacheck;
}

=head2 download_data

Downloads the data given a url and a destination for the data

=cut

sub download_data {
	my %args = @_;
	my $url         = $args{url};
	my $destination = $args{destination};
	`mkdir -p $destination`;

	# FIXME this is insecure!
	# but then again, so is just about every other line of code in this program...
	my $ff = File::Fetch->new( uri => $url );
	$ff->fetch( to => "$destination/.." );
	debug "URL : $url\n";
	$url =~ /\/(\w+)\.tgz/;
	my $archive = $1;
	debug "ARCHIVE : $archive\n";
	if ( -e "$destination/.. " ) {
		`rm -rf $destination/..`;
	}
	`cd $destination/../ ; tar xzf $archive.tgz ; touch $archive`;
	`rm $destination/../$archive.tgz`;
}
my $marker_update_url           = "http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift_markers/markers.tgz";
my $markers_extended_update_url = "http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift_markers/markers_extended.tgz";
my $ncbi_url                    = "http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz";

=head2 marker_update_check

Checks with the remote location to see if the markers have been updated since the last download
If the downloadable marker package is newer than what is currently being used then download the newer version.

=cut

sub marker_update_check {
	my %args        = @_;
	my $self        = $args{self};
	my $url         = $args{url};
	my $marker_path = $args{dir};
	my ( $content_type, $document_length, $modified_time, $expires, $server ) = head($url);
	debug "MARKER_PATH : " . $marker_path . "\n";
	my $get_new_markers = 0;
	if ( -x $marker_path ) {
		my $mtime = ( stat($marker_path) )[9];
		debug "TEST LOCAL :" . localtime($mtime) . "\n";
		if ( !defined($modified_time) ) {
			warn "Warning: unable to connect to marker update server, please check your internet connection\n";
		} elsif ( $modified_time > $mtime ) {
			debug "TEST REMOTE:" . localtime($modified_time) . "\n";
			warn "Found newer version of the marker data\n";
			$get_new_markers = 1;
		}
	} else {
		if ( !defined($modified_time) ) {
			croak "Marker data not found and unable to connect to marker update server, please check your phylosift configuration and internet connection!\n";
		}
		warn "Unable to find marker data!\n";
		$get_new_markers = 1;
	}
	if ($get_new_markers) {
		warn "Downloading from $url\n";
		download_data( url=>$url, destination=>$marker_path );
		my @markers = gather_markers( self => $self, path => $marker_path );
		index_marker_db( self => $self, markers => \@markers, path => $marker_path );
	}
}

=head2 data_checks

Performs checks to see if the NCBI data available for download is newer than what is currently being used by PhyloSift
If the available data is newer than the used data, download the newer version.

=cut

sub data_checks {
	my %args = @_;
	my $self = $args{self};

	#
	# first check for the standard marker directory
	$marker_dir = get_data_path( data_name=>"markers", data_path=>$Phylosift::Settings::marker_path );
	marker_update_check( self => $self, dir => $marker_dir, url => $marker_update_url );

	#
	# now look for the extended marker directory
	$markers_extended_dir = get_data_path(data_name=> "markers_extended", data_path=>$Phylosift::Settings::marker_path );
	if ( $self->{"extended"} ) {
		marker_update_check( self => $self, dir => $markers_extended_dir, url => $markers_extended_update_url );
	}

	#
	# now check for the NCBI taxonomy data
	$ncbi_dir = get_data_path(data_name=> "ncbi", data_path=>$Phylosift::Settings::ncbi_path );
	my ( $content_type, $document_length, $modified_time, $expires, $server ) = head("$ncbi_url");
	if ( -x $ncbi_dir ) {
		my $ncbi_time = ( stat($ncbi_dir) )[9];
		if ( !defined($modified_time) ) {
			warn "Warning: unable to connect to NCBI taxonomy update server, please check your internet connection\n";
		} elsif ( $modified_time > $ncbi_time ) {
			warn "Found newer version of NCBI taxonomy data!\n";
			warn "Downloading from $ncbi_url\n";
			download_data( url=>$ncbi_url, destination=>$ncbi_dir );
		}
	} else {
		if ( !defined($modified_time) ) {
			croak "NCBI taxonomy data not found and unable to connect to update server, please check your phylosift configuration and internet connection!\n";
		}
		warn "Unable to find NCBI taxonomy data!\n";
		warn "Downloading from $ncbi_url\n";
		download_data( url=>$ncbi_url, destination=>$ncbi_dir );
	}
}

=head2 fasta2stockholm

Convert a bunch of fasta files to stockholm format
This code is adapted from fasta2stockholm.pl, (c) Ian Holmes and licensed under the GPL
See https://github.com/ihh/dart/ for the original source code

=cut

sub fasta2stockholm {
	my %args = @_;
	my $fasta  = $args{fasta};
	my $output = $args{output};
	open( STOCKOUT, ">$output" );

	# read FASTA file
	my @seq;
	my @name;
	my $name;
	my $curseq = "";
	open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
	while (<FASTA>) {
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
	close FASTA;

	# check all seqs are same length
	my $length;
	my $lname;
	for ( my $sI = 0 ; $sI < @name ; $sI++ ) {
		my $nname = $name[$sI];
		my $sseq  = $seq[$sI];
		my $l     = length $sseq;
		if ( defined $length ) {
			croak "Sequences not all same length ($lname is $length, $nname is $l)" unless $length == $l;
		} else {
			$length = length $sseq;
			$lname  = $nname;
		}
	}

	# print Stockholm output
	print STOCKOUT "# STOCKHOLM 1.0\n";
	for ( my $sI = 0 ; $sI < @name ; $sI++ ) {
		print STOCKOUT $name[$sI], " ", $seq[$sI], "\n";
	}
	print STOCKOUT "//\n";
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
	my $seq_line = 0;	# counter for which line we're on in each block.
	while (my $line = <$STREAMIN>) {
		if( $line  !~ /\S/ || $line =~ /^\s*#/ || length($line)<2){
			$seq_line = 0;
			next;
		}
		if ($line =~ /^\s*\/\//) { $outbuffer .= printseq( columns => $columns, sorted => $sorted, seq => \@seq, names => \@names, out => $args{out} ); }
		else {
			chomp $line;
			my ( $name, $seq ) = split /\s+/, $line;
			$seq =~ s/[\.\-]//g unless $gapped;
			$seq[$seq_line] .= $seq;
			$names[$seq_line] = $name;
			$seq_line++;
		}
	}
	$outbuffer .= printseq( columns => $columns, sorted => $sorted, seq => \@seq, names => \@names, out => $args{out} );
	return $outbuffer;
}

sub printseq {
	my %args = @_;
	my @seq = @{$args{seq}};
	my @names = @{$args{names}};
	my $out  = "";
	for(my $j=0; $j<@names; $j++){
		$out .= ">$names[$j]\n";
		for ( my $i = 0 ; $i < length( $seq[$j] ) ; $i += $args{columns} ) {
			$out .= substr( $seq[$j], $i, $args{columns} ) . "\n";
		}
	}
	return $out;
}

=head2 read_name_table

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub read_name_table {
	my %args=@_;
	my $markerDir = $args{marker_directory};
	my %result;
	open( ALINAMES, "grep \">\" $markerDir/*.ali |" );
	while ( my $line = <ALINAMES> ) {
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
	my %args = @_;
	my $self   = $args{self};
	my $marker = $args{marker};
	my $length = 0;
	if ( is_protein_marker( marker => $marker ) ) {
		my $hmm_file = get_marker_hmm_file( self=>$self, marker=>$marker );
		open( HMM, $hmm_file ) || croak "Unable to open $hmm_file\n";
		while ( my $line = <HMM> ) {
			if ( $line =~ /LENG\s+(\d+)/ ) {
				$length = $1;
				last;
			}
		}
	} else {
		my $cm_file = get_marker_cm_file( self=>$self, marker=>$marker );
		open( CM, $cm_file ) || croak "Unable to open $cm_file\n";
		while ( my $line = <CM> ) {
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
	my %args = @_;
	my $self          = $args{self};
	my $marker        = $args{marker};
	my $gapmultiplier = $args{gap_multiplier};
	my $len           = get_marker_length( $self, $marker );
	my $glen;
	$glen = "P" x $len x $gapmultiplier if $gapmultiplier == 1;
	$glen = "A" x $len x $gapmultiplier if $gapmultiplier == 3;
	my $newseq = Bio::LocatableSeq->new( -seq => $glen, -id => "dummydummydummy", start => 1, end => ( $len * $gapmultiplier ) );
	my $aln = Bio::SimpleAlign->new();
	$aln->add_seq($newseq);
	return $aln;
}

=head2 get_marker_path

Determines the filesystem path to a marker. Searches for the marker in the base marker directory, the extended markers, and any local markers

=cut

sub get_marker_path {
	my %args   = @_;
	my $self   = $args{self};
	my $marker = $args{marker};
	
	# check for old-style marker first
	return "$marker_dir" if ( -e "$marker_dir/$marker.faa" );

	# check for new-style in standard directory
	return "$marker_dir" if ( -d "$marker_dir/$marker" );

	# check for new-style in extended directory
	return "$markers_extended_dir" if ( -d "$markers_extended_dir/$marker" );

	# TODO: check any local marker repositories
	warn "Could not find repository for marker $marker\n";
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

=head2 get_alignment_marker_file 

Returns the alignment file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_alignment_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.ali";
	} else {
		return "$marker_path/$bname.updated.fasta";
	}
}

=head2 get_marker_aln_file

Returns the aligned fasta file for the marker

=cut

sub get_marker_aln_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.ali" if ( -e "$marker_path/$bname.ali" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.aln";
	} else {
		return "$marker_path/$marker.updated/$bname.ali";
	}
}

=head2 get_marker_rep_file

Returns the fasta file of unaligned full length representative sequences for the marker

=cut

sub get_marker_rep_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.faa" if ( -e "$marker_path/$bname.faa" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.rep";
	} else {
		return "$marker_path/$marker.updated/$bname.rep";
	}
}

=head2 get_marker_hmm_file

Returns the HMM file for the marker

=cut

sub get_marker_hmm_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $local       = $args{loc} || 0; #using loc instead of local (reserved word)
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.hmm" if ( -e "$marker_path/$bname.hmm" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.hmm";
	} else {
		return "$marker_path/$bname.hmm";
	}
}

=head2 get_marker_cm_file

Returns the CM (infernal covarion model) file for the marker

=cut

sub get_marker_cm_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	my $updated     = $self->{"updated"} ? ".updated" : "";
	return "$marker_path/$marker$updated/$bname.cm";
}

=head2 get_marker_stockholm_file

Returns the stockholm file for the marker

=cut

sub get_marker_stockholm_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.stk" if ( -e "$marker_path/$bname.ali" );

		# using new-style marker directories
		return "$marker_path/$marker/$bname.stk";
	} else {
		return "$marker_path/$marker.updated/$bname.stk";
	}
}

=head2 get_marker_taxon_map

Returns the path to the lookup table between marker gene IDs and their taxa

=cut

sub get_marker_taxon_map {
	my %args = @_;
	my $self = $args{self};
	return "$marker_dir/marker_taxon_map.updated.txt" if ( $self->{"updated"} );
	return "$marker_dir/marker_taxon_map.txt";
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
	return 1 if ( -e "$marker_path/$marker/$bname.hmm" );
	return 1 if ( -e "$marker_path/$bname.hmm" );
	return 0 if ( -e "$marker_path/$marker/$bname.cm" );
	return 1;
}

=head2 get_fasta_marker_file

Returns the fasta file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_fasta_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.faa";
	} else {
		return "$marker_path/$bname.faa";
	}
}

=head2 get_marker_package

Returns the path to the marker package

=cut

sub get_marker_package {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$marker";
	} else {
		return "$marker_path/$marker.updated";
	}
}

=head2 get_aligner_output_fasta_AA
Returns the FastA file containing amino acid read or contig alignments to the marker
given by markerName
=cut

sub get_aligner_output_fasta_AA {
	my %args = @_;
	my $marker = $args{marker};
	my $bname = get_marker_basename( marker => $marker );
	return "$bname.trim.fasta";
}

=head2 get_aligner_output_fasta_DNA
Returns the FastA file containing DNA read or contig alignments to the marker
given by markerName
=cut

sub get_aligner_output_fasta_DNA {
	my %args = @_;
	my $marker = $args{marker};
	my $bname = get_marker_basename( marker => $marker );
	return "$bname.trim.fna.fasta";
}

=head2 get_read_placement_file
Returns the read placement Jplace file to the marker
given by markerName
=cut

sub get_read_placement_file {
	my %args = @_;
	my $marker = $args{marker};
	my $bname = get_marker_basename( marker => $marker );
	return "$bname.trim.jplace";
}

=head2 get_read_placement_file_DNA
Returns the read placement Jplace file to the marker
given by markerName
=cut
sub get_read_placement_file_DNA {
	my %args = @_;
	my $marker = $args{marker};
	my $bname = get_marker_basename( marker => $marker );
	return "$bname.trim.fna.jplace";
}

=head2 get_trimfinal_marker_file

Returns the .trimfinal file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub get_trimfinal_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.trimfinal" if -e "$marker_path/$bname.trimfinal";
		return "$marker_path/$marker/$bname.aln";
	} else {
		return "$marker_path/$bname.trimfinal";
	}
}

=head2 get_trimfinal_fasta_marker_file

Returns the .trimfinal.fasta file for the markerName passed in as an argument
If the user chooses the updated markers, the updated file is returned instead

=cut

sub get_trimfinal_fasta_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.trimfinal.fasta" if -e "$marker_path/$bname.trimfinal.fasta";
		return "$marker_path/$marker/$bname.aln";
	} else {
		return "$marker_path/$bname.trimfinal.fasta";
	}
}

=head2 get_tree_marker_file

Returns the .final.tre file from the marker directory 
The user chooses the updated or stock version

=cut

sub get_tree_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.final.tre" if -e "$marker_path/$bname.final.tre";
		return "$marker_path/$marker/$bname.tre";
	} else {
		return "$marker_path/$bname.updated.tre";
	}
}

=head2 get_tree_stats_marker_file

Return the updated or stock version of the Tree stats file

=cut

sub get_tree_stats_marker_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
		return "$marker_path/$bname.in_phyml_stats.txt";
	} else {
		return "$marker_path/$bname.updated.RAxML_info";
	}
}

=head2 get_ncbi_map_file

Returns the updated of stock version of the NCBI map file

=cut

sub get_ncbi_map_file {
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_path = get_marker_path( self => $self, marker => $marker );
	my $bname       = get_marker_basename( marker => $marker );
	if ( $self->{"updated"} == 0 ) {
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
	my %args = @_;
	my $self        = $args{self};
	my $marker      = $args{marker};
	my $marker_file = get_marker_rep_file( self=>$self, marker=>$marker );
	my $rep_num     = `grep -c '>' $marker_file`;
	chomp($rep_num);
	return $rep_num;
}

=head2 concatenate_alignments

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub concatenate_alignments {
	my %args = @_;
	my $self          = $args{self};
	my $outputFasta   = $args{output_fasta};
	my $outputMrBayes = $args{output_bayes};
	my $gapmultiplier = $args{gap_multiplier};    # 1 for protein, 3 for reverse-translated DNA
	my @alignments    = $args{alignments};
	my $catobj        = 0;
	open( MRBAYES, ">$outputMrBayes" );
	my $partlist = "partition genes = " . scalar(@alignments) . ": ";
	my $prevlen  = 0;

	foreach my $file (@alignments) {
		my $aln;
		my $marker = basename($file);
		$marker =~ s/\..+//g;     # FIXME: this should really come from a list of markers
		unless ( -e $file ) {

			# this marker doesn't exist, need to create a dummy with the right number of gap columns
			$aln = make_dummy_file(self=> $self, marker=>$marker,gap_multiplier=> $gapmultiplier );
		} else {
			my $in = Bio::AlignIO->new( -file => $file, '-format' => 'fasta' );
			unless ( $aln = $in->next_aln() ) {

				# empty marker alignment file, need to create a dummy with the right number of gap columns
				$aln = make_dummy_file( self=>$self, marker=>$marker, gap_multiplier=>$gapmultiplier );
			}
		}
		my $csname = $file;
		$csname =~ s/\..+//g;
		$partlist .= "," if $catobj != 0;
		$partlist .= $csname;
		print MRBAYES "charset $csname = " . ( $prevlen + 1 ) . "-" . ( $prevlen + $aln->length() ) . "\n";
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
				my $newseq = Bio::LocatableSeq->new( -seq => $tmpseq, -alphabet => "protein", -id => $catseq->id, start => 0, end => 0 );
				$newseq->verbose(-1);
				$aln->add_seq($newseq);
			}
		}

		# vice versa
		foreach my $alnseq ( $aln->each_alphabetically() ) {
			if ( length( $catobj->each_seq_with_id( $alnseq->id ) == 0 ) ) {

				# add this sequence as all gaps
				my $tmpseq = "-" x $catobj->length();
				my $newseq = Bio::LocatableSeq->new( -seq => $tmpseq, -alphabet => "protein", -id => $alnseq->id, start => 0, end => 0 );
				$newseq->verbose(-1);
				$catobj->add_seq($newseq);
			}
		}
		$catobj = cat( $catobj, $aln );
		$catobj->verbose(-1);
	}
	return unless $catobj;	# exit if there's nothing to write about
	
	print MRBAYES "$partlist;\n";
	my $out = Bio::AlignIO->new( -file => ">$outputFasta", '-format' => 'fasta' );
	foreach my $dummyseq ( $catobj->each_seq_with_id("dummydummydummy") ) {
		$catobj->remove_seq($dummyseq);
	}
	$out->write_aln($catobj);
}
my %timers;

sub start_timer {
	my %args= @_;
	my $timername = $args{name};
	my $t         = time;
	my @timerval  = localtime($t);
	$timers{$timername} = $t;
	debug sprintf( "Before $timername %4d-%02d-%02d %02d:%02d:%02d\n",
				   $timerval[5] + 1900,
				   $timerval[4] + 1,
				   $timerval[3], $timerval[2], $timerval[1], $timerval[0] );
}

sub end_timer {
	my %args = @_;
	my $timername = $args{name};
	my $t         = time;
	my @timerval  = localtime($t);
	debug sprintf( "After $timername %4d-%02d-%02d %02d:%02d:%02d\n",
				   $timerval[5] + 1900,
				   $timerval[4] + 1,
				   $timerval[3], $timerval[2], $timerval[1], $timerval[0] );
	return $t - $timers{$timername};
}

=head2 marker_oldstyle

Checks whether a marker is in the old style format (PMPROK*) or the new format (marker package directory)

=cut

sub marker_oldstyle {
	my %args = @_;
	my $marker = $args{markers};
	return 1 if ( $marker =~ /PMPROK/ );
	return 0;
}

=head2 open_SeqIO_object

Opens a sequence file and returns a SeqIO object.  Allows for gzip and bzip compression
returns a SeqIO object

=cut

sub open_SeqIO_object {
	my %args = @_;
	my $format = $args{format} || "FASTA";    #default
	my $io_object;
	if ( exists $args{format} ) {
		$format = $args{format};
	}
	if ( $args{file} =~ /\.gz$/ ) {
		$io_object = Bio::SeqIO->new( -file => "zcat $args{file} |", -format => $format );
	} elsif ( $args{file} =~ /\.bz2$/ ) {
		$io_object = Bio::SeqIO->new( -file => "bzcat $args{file} |", -format => $format );
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
	my $F1IN;
	if ( $args{file} =~ /\.gz$/ ) {
		open( $F1IN, "zcat $args{file} |" ) or croak "Can't open $args{file}\n";
	} elsif ( $args{file} =~ /\.bz2$/ ) {
		open( $F1IN, "bzcat $args{file} |" );
	} else {
		open( $F1IN, $args{file} );
	}
	return $F1IN;
}

sub get_db {
	my %args    = @_;
	my $path    = $args{path} || $marker_dir;
	my $self    = $args{self};                  # optional argument;
	my $db_name = $args{db_name};
	if ( defined($self) && !defined( $args{path} ) ) {
		return $markers_extended_dir . "/$db_name" if defined( $self->{"extended"} && $self->{"extended"} );
		return $marker_dir . "/$db_name";
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

sub get_rapsearch_db {
	my %args = @_;
	$args{db_name} = "rep";
	return get_db(%args);
}

=head2 get_bowtie2_db

returns the name and path of the bowtie2 DB

=cut

sub get_bowtie2_db {
	my %args = @_;
	$args{db_name} = "rnadb";
	return get_db(%args);
}

sub get_candidate_file {
	my %args   = @_;
	my $self   = $args{self};
	my $marker = $args{marker};
	my $type   = $args{type};
	my $dna    = $args{dna};
	my $new    = $args{new};
	my $ffn    = "";
	$ffn = ".ffn" if ( defined( $args{dna} ) && $dna );
	my $candidate = ".candidate";
	$candidate = ".newCandidate" if defined( $args{new} ) && $new;
	my $dir = $self->{"blastDir"};
	$dir = $self->{"alignDir"} if defined( $args{new} ) && $new;
	$marker =~ s/.+\///g;    # strip off any prepended directories
	return "$dir/$marker$type$candidate$ffn";
}

=head2 index_marker_db

Indexes the marker database for searches with rapsearch2, blastall, and bowtie2
Input: marker list and self

=cut

sub index_marker_db {
	my %args    = @_;
	my $self    = $args{self};
	my $path    = $args{path};
	my @markers = @{ $args{markers} };

	# use alignments to make an unaligned fasta database containing everything
	# strip gaps from the alignments
	my $bowtie2_db = get_bowtie2_db( path => $path );
	my $bowtie2_db_fasta = "$bowtie2_db.fasta";
	open( my $PDBOUT, ">" . get_blastp_db( path => $path ) );
	open( my $RNADBOUT, ">" . $bowtie2_db_fasta );
	foreach my $marker (@markers) {
		my $marker_rep = get_marker_rep_file( self=>$args{self}, marker=>$marker );
		my $DBOUT = $RNADBOUT;
		$DBOUT = $PDBOUT if is_protein_marker( marker => $marker );
		unless ( -f $marker_rep ) {
			warn "Warning: marker $marker appears to be missing data\n";
			next;
		}
		open( INALN, $marker_rep );
		while ( my $line = <INALN> ) {
			if ( $line =~ /^>(.+)/ ) {
				print $DBOUT "\n>$marker" . "__$1\n";
			} else {
				$line =~ s/[-\.\n\r]//g;
				$line =~ tr/a-z/A-Z/;
				print $DBOUT $line;
			}
		}
	}
	print $PDBOUT "\n";
	print $RNADBOUT "\n";
	close $PDBOUT;    # be sure to flush I/O
	close $RNADBOUT;
	#print $path."\n";
	# make a blast database
	my $blastp_db = get_blastp_db( path => $path );
	system("$Phylosift::Utilities::formatdb -i $blastp_db -o F -p T -t RepDB");
	`mv $blastp_db $path/rep.dbfasta`;

	# make a rapsearch database
	`cd $path ; $Phylosift::Utilities::preRapSearch -d rep.dbfasta -n rep`;

	# make a last database
	`cd $path ; $Phylosift::Utilities::lastdb -p -c replast rep.dbfasta`;
	unlink("$path/rep.dbfasta");    # don't need this anymore!

	# make a bowtie2 database
	if ( -e "$bowtie2_db_fasta" ) {
		`cd $path ; $Phylosift::Utilities::lastdb -c $bowtie2_db $bowtie2_db_fasta`;
		`cd $path ; $Phylosift::Utilities::bowtie2build $bowtie2_db_fasta $bowtie2_db`;
	}
	
	# now create the .hmm files if they aren't already present
	# this is the case in the extended marker set, since the hmms are too big for transit
	foreach my $marker (@markers) {
		my $hmm_file = get_marker_hmm_file($args{self}, $marker);
		my $cm_file = get_marker_cm_file($args{self}, $marker);
		next if -e $hmm_file;
		next if -e $cm_file;
		next if is_protein_marker( marker => $marker );
		my $stk_file = get_marker_stockholm_file($args{self}, $marker);
		`$hmmbuild $hmm_file $stk_file`;
	}
}

=head2 gather_markers

=item *

    ARGS : $markerFile - file to read the marker names from that will be used in the pipeline.
    Reads markers from a file if specified by the user OR gathers all the markers from the markers directory
    The Marker names are stored in an Array
    If the filename is empty, use all the default markers.

=back

=cut

sub gather_markers {
	my %args        = @_;
	my $self        = $args{self};
	my $marker_file = $args{marker_file};
	my $path        = $args{path} || $marker_dir;
	my @marks       = ();

	#create a file with a list of markers called markers.list
	if ( defined($marker_file) && -f $marker_file && $marker_file ne "" ) {

		#gather a custom list of makers, list convention is 1 marker per line
		open( markersIN, $marker_file );
		while (<markersIN>) {
			chomp($_);
			push( @marks, $_ );
		}
		close(markersIN);
	} else {

		# gather all markers
		# this is for the original marker set
		my @files = <$path/*.faa>;
		foreach my $file (@files) {
			$file =~ m/\/(\w+).faa/;
			push( @marks, $1 );
		}

		# now gather directory packaged markers (new style)
		# use maxdepth 2 for two-directory-level markers in the extended marker set
		open( MLIST, "find $path -maxdepth 2 -mindepth 1 -type d |" );
		while ( my $line = <MLIST> ) {
			chomp $line;
			next if $line =~ /PMPROK/;
			next if $line =~ /concat/;
			next if $line =~ /representatives/;
			$line = substr( $line, length($path) + 1 );

			# all markers need to have an hmm or a cm else they are not usable
			my $baseline = $line;
			$baseline =~ s/.+\///g;
			next unless ( -e "$path/$line/$line.cm" || -e "$path/$line/$baseline.hmm" );
			push( @marks, $line );
		}
	}
	return @marks;
}

=head2 get_sequence_input_type

Checks whether input is FastA, FastQ, which quality type (33 or 64), and DNA or AA
Returns a hash reference with the values 'seqtype', 'format', and 'qtype' populated.

=cut

sub get_sequence_input_type {
	my $file = shift;
	my %type;
	my $FILE       = open_sequence_file( file => $file );
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

	while ( my $line = <$FILE> ) {
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
			$counter  += length($line) - 1;
			$dnacount += $line =~ tr/[ACGTNacgtn]//;
			$allcount += length($line) - 1;
		}
		$line_count++;
		last if ( $line_count > 1000 );
	}
	close($FILE);
	$maxfound = $counter > $maxfound ? $counter : $maxfound;
	$type{seqtype} = "protein" if ( $dnacount < $allcount * 0.75 );
	$type{seqtype} = "dna"     if ( $type{format} eq "fastq" );       # nobody using protein fastq (yet)
	$type{qtype}   = "phred64" if $minq < 255;
	$type{qtype}   = "phred33" if $minq < 64;
	$type{paired} = 0;                                                # TODO: detect interleaved read pairing
	return \%type;
}

=head2 get_sequence_input_type_quickndirty

Checks whether input is either short sequence reads, e.g. < 500nt or assembled fragments
without reading the whole file.

=cut

sub get_sequence_input_type_quickndirty {
	my $file         = shift;
	my $maxshortread = 500;
	open( FILE, $file );
	my $filesize = -s "$file";
	my $counter  = 0;
	my $maxfound = 0;
	my $allcount = 0;
	my $dnacount = 0;
	my $seqtype  = "dna";
	my $length   = "long";
	my $format   = "unknown";
	for ( my $i = 0 ; $i < 200 ; $i++ ) {
		my $seekpos = int( rand( $filesize - 100 ) );
		$seekpos = 0 if ( $i == 0 );    # always start with the first line in case the sequence is on a single line!
		seek( FILE, $seekpos, 0 );
		$counter = 0;
		my $line = <FILE>;              # burn a line to be sure we get to sequence
		while ( $line = <FILE> ) {
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
	$seqtype = "dna" if ( $format eq "fastq" );    # nobody using protein fastq (yet)
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
	$datestr .= 0 if $timerval[4] <= 9;
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
	print "(c) 2011, 2012 Aaron Darling and Guillaume Jospin\n";
	print "\nCITATION:\n";
	print "		PhyloSift. A. Darling, H. Bik, G. Jospin, J. A. Eisen. Manuscript in preparation\n";
	print "\n\nPhyloSift incorporates several other software packages, please consider also citing the following papers:\n";
	print qq{

		pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree.
		Frederick A Matsen, Robin B Kodner, and E Virginia Armbrust
		BMC Bioinformatics 2010, 11:538
		
		RAPSearch2: a fast and memory-efficient protein similarity search tool for next generation sequencing data.
		Yongan Zhao, Haixu Tang, and Yuzhen Ye
		Bioinformatics (2011)
		
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
	my %args = @_;
	my $aln_file   = $args{aln};
	my $target_dir = $args{target_dir};
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $in = open_SeqIO_object( file => $aln_file );
	open( FILEOUT, ">" . $target_dir . "/" . $core . ".fasta" ) or carp("Couldn't open $target_dir$core.fasta for writing\n");
	while ( my $seq_object = $in->next_seq() ) {
		my $seq = $seq_object->seq;
		my $id  = $seq_object->id;
		$seq =~ s/-\.//g;    # shouldnt be any gaps
		print FILEOUT ">" . $id . "\n" . $seq . "\n";
	}
	close(FILEOUT);
	return "$target_dir/$core.fasta";
}

=head2 generate_hmm

input: alignment_file, target_directory
generates a HMM profile from an alignement in FASTA format (arg) using hmmbuild.  The hmm is placed in the target_directory

=cut

sub generate_hmm {
	my %args=@_;
	my $file_name  = $args{file};
	my $target_dir = $args{target_dir};
	my ( $core_name, $path, $ext ) = fileparse( $file_name, qr/\.[^.]*$/ );
	if ( !-e "$target_dir/$core_name.hmm" ) {
		`hmmbuild --informat afa $target_dir/$core_name.hmm $file_name`;
	}
	return "$target_dir/$core_name.hmm";
}

=head2 hmmalign_to_model

input : hmm_profile,sequence_file,target_dir
Aligns sequences to an HMM model and outputs an alignment
=cut

sub hmmalign_to_model {
	my %args=@_;
	my $hmm_profile   = $args{hmm};
	my $sequence_file = $args{sequence_file};
	my $target_dir    = $args{target_dir};
	my $ref_ali       = $args{reference_aln};
	my ( $core_name, $path, $ext ) = fileparse( $sequence_file, qr/\.[^.]*$/ );
	if ( !-e "$target_dir/$core_name.aln" ) {
		`hmmalign --mapali $ref_ali --trim --outformat afa -o $target_dir/$core_name.aln $hmm_profile $sequence_file`;

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
		my %args = @_;
	my $aln_file    = $args{aln};
	my $output_path = $args{output_path};
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $in = open_SeqIO_object( file => $aln_file );
	my $seq_count = 0;
	open( FILEOUT, ">$output_path" ) or carp("Couldn't open $output_path for writing\n");
	while ( my $seq_object = $in->next_seq() ) {
		my $seq = $seq_object->seq;
		my $id  = $seq_object->id;
		$seq =~ s/[-\.]//g;    # shouldnt be any gaps
		print FILEOUT ">" . $id . "\n" . $seq . "\n";
		$seq_count++;
	}
	close(FILEOUT);
	return $seq_count;
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
