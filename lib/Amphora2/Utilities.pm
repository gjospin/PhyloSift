package Amphora2::Utilities;

use 5.006;
use strict;
use warnings;
use Carp;
use File::Basename;
use Bio::AlignIO;
use Bio::Align::Utilities qw(:all); 
use POSIX ();
use LWP::Simple;
use File::Fetch;
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
our $rapSearch= "";
our $preRapSearch = "";
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

	$rapSearch = get_program_path("rapsearch",$Amphora2::Settings::a2_path);
	if($rapSearch eq ""){
	    carp("rapsearch was not found\n");
	    return 1;
	}

	$preRapSearch = get_program_path("prerapsearch",$Amphora2::Settings::a2_path);


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
	if( defined($datapath) && $datapath ne "" ){
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
	my $ff = File::Fetch->new(uri=>$url);
	$ff->fetch(to=>"$destination/..");
#	`wget -N -nv $url -O  $destination/../amphora_data.tar.gz`;
	print "URL : $url\n";
	$url =~ /\/(\w+)\.tgz/;
	my $archive = $1;
	print "ARCHIVE : $archive\n";
	if(-e "$destination/.. "){
	    `rm -rf $destination/..`;
	}
	`cd $destination/../ ; tar xzf $archive.tgz ; touch $archive`;
	`rm $destination/../$archive.tgz`;
}


my $marker_update_url = "http://edhar.genomecenter.ucdavis.edu/~mlangill/markers.tgz";
my $ncbi_url = "http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz";

sub dataChecks {
	$marker_dir = get_data_path( "markers", $Amphora2::Settings::marker_path );
	my ($content_type, $document_length, $modified_time, $expires, $server)= head("$marker_update_url");
	print "TEST REMOTE:".localtime($modified_time)."\n";
	print STDERR "MARKER_PATH : ".$marker_dir."\n";
	if(-x $marker_dir){
	    my $mtime = (stat($marker_dir))[9];
	    print "TEST LOCAL :".localtime($mtime)."\n";	
	    if($modified_time > $mtime){
		warn "Found newer version of the marker data\n";
		warn "Downloading from $marker_update_url\n";
		download_data( $marker_update_url, $marker_dir );
	    }
	}else{
	    warn "Unable to find marker data!\n";
	    warn "Downloading from $marker_update_url\n";
	    download_data($marker_update_url, $marker_dir);
	}
	$ncbi_dir = get_data_path( "ncbi", $Amphora2::Settings::ncbi_path );
	($content_type, $document_length, $modified_time, $expires, $server)= head("$ncbi_url");
	if( -x $ncbi_dir ){
	    my $ncbi_time =(stat($ENV{"HOME"}."/share/amphora2/ncbi"))[9];
	    if($modified_time > $ncbi_time){
		warn "Found newer version of NCBI taxonomy data!\n";
		warn "Downloading from $ncbi_url\n";
		download_data( $ncbi_url, $ncbi_dir );
	    }
	}else{
	    warn "Unable to find NCBI taxonomy data!\n";
	    warn "Downloading from $ncbi_url\n";
	    download_data( $ncbi_url, $ncbi_dir);
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

=head2 concatenateAlignments

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub concatenateAlignments {
	my $outputFasta = shift;
	my $outputMrBayes = shift;
	my @alignments = @_;
	my $catobj = 0;
	open(MRBAYES, ">$outputMrBayes");
	my $partlist = "partition genes = ".scalar(@alignments).": ";
	my $prevlen = 0;
	foreach my $file(@alignments){
	    print STDERR $file."\n";
		my $in  = Bio::AlignIO->new(-file => $file , '-format' => 'fasta');
		while ( my $aln = $in->next_aln() ) {
			my $csname = $file;
			$csname =~ s/\..+//g;
			$partlist .= "," if $catobj != 0;
			$partlist .= $csname;
			print MRBAYES "charset $csname = ".($prevlen+1)."-".($prevlen+$aln->length())."\n";
			$prevlen += $aln->length();
			my $prevseq = 0;
			my $newaln = $aln->select(1,1);
			foreach my $curseq( $aln->each_alphabetically() ){
				if($prevseq==0){
					$prevseq = $curseq;
					next;
				}
				if($prevseq->id ne $curseq->id){
					$newaln->add_seq($curseq);
				}
				$prevseq = $curseq;
			}
			$aln = $newaln;

			if($catobj == 0){
				$catobj = $aln;
				next;
			}
			# add any sequences missing from this sequence
			foreach my $catseq ( $catobj->each_alphabetically() ){
				if(length( $aln->each_seq_with_id( $catseq->id ) == 0) ){
					# add this sequence as all gaps
					my $tmpseq = "-" x $aln->length();
					my $newseq = Bio::LocatableSeq->new( -seq => $tmpseq, -id => $catseq->id, start=>0, end=>$aln->length());
					$aln->add_seq($newseq);
			}
			}
			# vice versa
			foreach my $alnseq ( $aln->each_alphabetically() ){
				if(length( $catobj->each_seq_with_id( $alnseq->id ) == 0) ){
					# add this sequence as all gaps
					my $tmpseq = "-" x $catobj->length();
					my $newseq = Bio::LocatableSeq->new( -seq => $tmpseq, -id => $alnseq->id, start=>0, end=>$catobj->length());
					$catobj->add_seq($newseq);
				}
			}
			$catobj = cat( $catobj, $aln );
		}
	}
	print MRBAYES "$partlist;\n";
	my $out = Bio::AlignIO->new(-file => ">$outputFasta" , '-format' => 'fasta');
	$out->write_aln($catobj);
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
