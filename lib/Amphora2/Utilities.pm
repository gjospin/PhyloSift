package Amphora2::Utilities;

use 5.006;
use strict;
use warnings;
use Carp;
use File::Basename;

=head1 NAME

Amphora2::Utilities - Implements miscellaneous accessory functions for Amphora2

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

=head2 programChecks

checks the program requirements for Amphora-2
writes to STDERR if a program is missing or the wrong version is installed
returns 1 or 0 depending on success of failure.

=cut

sub get_program_path {
	my $progname = shift;
	my $progpath = shift;
	my $progcheck = "";
	if( defined($progpath) && $progpath ne "" && -x $progpath."/".$progname ){
		$progcheck = $progpath."/".$progname;
	}else{
		$progcheck = `which $progname`;
		chomp $progcheck;
	}
	return $progcheck;
}

# external programs used by Amphora2
our $pplacer = "";
our $hmmalign = "";
our $hmmsearch = "";
our $hmmbuild = "";
our $blastp = "";
our $makeblastdb = "";
our $translateSixFrame = "";
our $printneighbor = "";
our $Rscript = "";

sub programChecks {
	eval 'require Bio::Seq;';
	if ($@) {
	    carp "Bioperl was NOT found\n";
	    return 1;
	}

	$pplacer = get_program_path("pplacer", $Amphora2::Settings::pplacer_path);
	if($pplacer eq ""){
	    #program not found return;
	    carp("pplacer v1.1.alpha00 not found");
	    return 1;
	}elsif(`$pplacer -help` !~ m/v1.1.alpha00/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp("Warning : a different version of pplacer was found. Amphora-2 was tested with pplacer v1.1.alpha00\n");
	}

	$hmmalign = get_program_path("hmmalign", $Amphora2::Settings::hmmer3_path);
	if($hmmalign eq ""){
	    #program not found return;
	    carp("HMMER3 not found");
	    return 1;
	}elsif(`$hmmalign -h` !~ m/HMMER 3.0rc1/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp "Warning : a different version of HMMER was found. Amphora-2 was tested with HMMER 3.0rc1\n";
	}
	$hmmsearch = get_program_path("hmmsearch", $Amphora2::Settings::hmmer3_path);
	$hmmbuild = get_program_path("hmmbuild", $Amphora2::Settings::hmmer3_path);


	$blastp = get_program_path("blastp", $Amphora2::Settings::blast_path);
	if($blastp eq ""){
	    #program not found return;
	    carp("Blast 2.2.24+ not found");
	    return 1;
	}elsif(`$blastp -help` !~ m/BLAST 2.2.24+/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp "Warning : a different version of Blast was found. Amphora-2 was tested with BLAST 2.2.24+\n";
	}else{
	    #program found and the correct version is installed
	}

	# use makeblastdb from the same location as blastall
	$makeblastdb = get_program_path("makeblastdb", $Amphora2::Settings::blast_path);
	if($makeblastdb eq ""){
	    carp("makeblastdb from Blast+ not found");
	    return 1;
	}

	$translateSixFrame = get_program_path("translateSixFrame", $Amphora2::Settings::a2_path);
	if($translateSixFrame eq ""){
	    carp("Amphora2 translateSixFrame program not found");
	    return 1;
	}

	$printneighbor = get_program_path("printneighbor.R", $Amphora2::Settings::a2_path);
	if($printneighbor eq ""){
	    carp("Amphora2 printneighbor program not found");
	    return 1;
	}

	$Rscript = get_program_path("Rscript", $Amphora2::Settings::R_path);

	#ensure we have a place to put any R packages that might need to be downloaded
	my $amphora_r_lib_path = $ENV{"HOME"}."/.amphora2_rlibs";
	`mkdir -p $amphora_r_lib_path`;
	unless( defined($ENV{"R_LIBS_USER"}) && $ENV{"R_LIBS_USER"}=~/amphora2/ ){
		$ENV{"R_LIBS_USER"} .= ":$amphora_r_lib_path";
	}

	return 0;
}


=head2 dataChecks

Check for requisite Amphora-2 marker datasets

=cut

our $marker_dir = "";
our $ncbi_dir = "";

sub get_data_path {
	my $dataname = shift;
	my $datapath = shift;
	my $datacheck = "";
	if( defined($datapath) && $datapath ne "" && -x $datapath."/".$dataname ){
		$datacheck = $datapath."/".$dataname;
	}else{
		my $scriptpath = dirname($0);
		$scriptpath =~ s/bin\/?$//g;
		$datacheck = $scriptpath."/share/amphora2/".$dataname;
		return $datacheck if ( -x $datacheck );
		# if the system data dir doesn't exist, default to the user's home
		$datacheck = $ENV{"HOME"}."/share/amphora2/".$dataname;
	}
	return $datacheck;
}

sub download_data {
	my $url = shift;
	my $destination = shift;
	`mkdir -p $destination`;
	# FIXME this is insecure!
	`wget $url -O $destination/../amphora_data.tar.gz`;
	`cd $destination/../ ; tar xzf amphora_data.tar.gz`;
	`rm $destination/../amphora_data.tar.gz`;
}

my $marker_update_url = "http://edhar.genomecenter.ucdavis.edu/~mlangille/markers.tgz";
my $ncbi_url = "http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz";

sub dataChecks {
	$marker_dir = get_data_path( "markers", $Amphora2::Settings::marker_path );
	unless( -x $marker_dir ){
		warn "Unable to find marker data!\n";
		warn "Downloading from $marker_update_url\n";
		download_data( $marker_update_url, $marker_dir );
	}
	$ncbi_dir = get_data_path( "ncbi", $Amphora2::Settings::ncbi_path );
	unless( -x $ncbi_dir ){
		warn "Unable to find NCBI taxonomy data!\n";
		warn "Downloading from $ncbi_url\n";
		download_data( $ncbi_url, $ncbi_dir );
	}
}

=head2 fasta2stockholm

Convert a bunch of fasta files to stockholm format

=cut

sub fasta2stockholm {
	my $fasta = shift;
	my $output = shift;
	open( STOCKOUT, ">$output" );
	# read FASTA file
	    my %seq;
	    my @name;
	    my $name;
	    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
	    while (<FASTA>) {
		if (/^\s*>\s*(\S+)/) {
		    $name = $1;
		    croak "Duplicate name: $name" if defined $seq{$name};
		    push @name, $name;
		} else {
		    if (/\S/ && !defined $name) {
			warn "Ignoring: $_";
		    } else {
			s/\s//g;
			$seq{$name} .= $_;
		    }
		}
	    }
	    close FASTA;

	# check all seqs are same length
	    my $length;
	    my $lname;
	    foreach my $nname (@name) {
		my $l = length $seq{$nname};
		if (defined $length) {
		    croak "Sequences not all same length ($lname is $length, $nname is $l)" unless $length == $l;
		} else {
		    $length = length $seq{$nname};
		    $lname = $nname;
		}
	    }

	# print Stockholm output
	    print STOCKOUT "# STOCKHOLM 1.0\n";
	    foreach my $nname (@name) {
		print STOCKOUT $nname, " ", $seq{$nname}, "\n";
	    }
	    print STOCKOUT "//\n";
}

=head2 makeNameTable

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub readNameTable {
	my $markerDir = shift;
	my %result;
	open(ALINAMES, "grep \">\" $markerDir/*.ali |");
	while( my $line = <ALINAMES> ){
		$line =~ s/.+\:\>//g;
		my $commonName;
		if($line =~ /\{(.+)\}.+\[.+\]/){
			$commonName = $1;
		}elsif($line =~ /\[(.+)\]$/){
			$commonName = $1;
		}
		my @fields = split(/\s+/, $line);
		$result{$fields[0]}=$commonName;
	}
	return %result;
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

    perldoc Amphora2::Utilities


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

1; # End of Amphora2::Utilities
