package Amphora2::Utilities;

#use 5.006;
use strict;
use warnings;
use FindBin qw($Bin);
BEGIN { unshift(@INC, "$FindBin::Bin/legacy/") if $[ < 5.010; }
use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Align::Utilities qw(:all); 
use POSIX ();
use LWP::Simple;
use Carp;
require File::Fetch;
if($^O=~/arwin/){
	use lib "$FindBin::Bin/osx/darwin-thread-multi-2level/";
}

use Exporter;
use vars            qw[ @EXPORT @EXPORT_OK %EXPORT_TAGS @ISA ];;

@ISA            = 'Exporter';
@EXPORT         = qw[start_timer end_timer debug];
@EXPORT_OK      = qw[];

%EXPORT_TAGS    = (
STD     => \@EXPORT,
all     => [ @EXPORT, @EXPORT_OK ],
);        

our $debuglevel = 3;
sub debug {
	my $msg = shift;
	my $msglevel = shift || 2;
	print $msg if $debuglevel >= $msglevel;
}

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
#	print STDERR "BIN : $Bin\n";
#	exit;
	if( defined($progpath) && $progpath ne "" && -x $progpath."/".$progname ){
		$progcheck = $progpath."/".$progname;
	}else{
		$progcheck = `which $progname`;
		chomp $progcheck;
	}
	# last ditch attempt, check the directories from where the script is running
	$progcheck = $Bin."/".$progname unless( $progcheck =~ /$progname/ || !(-x $Bin."/".$progname) );
	$progcheck = $Bin."/bin/".$progname unless( $progcheck =~ /$progname/  || !(-x $Bin."/bin/".$progname) );
	# check the OS and use Mac binaries if needed
	if($^O =~ /arwin/){
		$progcheck = $Bin."/osx/".$progname unless( $progcheck =~ /$progname/  && !(-x $Bin."/".$progname) );
	}else{
	    #I don't think this line is needed because it is taken care of 5 lines before
		#$progcheck = $Bin."/bin/".$progname unless( $progcheck =~ /$progname/  && !(-x $Bin."/".$progname) );
	}
	return $progcheck;
}

# external programs used by Amphora2
our $pplacer = "";
our $hmmalign = "";
our $hmmsearch = "";
our $hmmbuild = "";
our $blastall = "";
our $formatdb = "";
our $rapSearch= "";
our $preRapSearch = "";
our $raxml = "";
sub programChecks {
	eval 'require Bio::Seq;';
	if ($@) {
	    carp "Bioperl was NOT found\n";
	    return 1;
	}

	$pplacer = get_program_path("pplacer", $Amphora2::Settings::pplacer_path);
	if($pplacer eq ""){
	    #program not found return;
	    carp("pplacer v1.1.alpha09 not found");
	    return 1;
	}elsif(`$pplacer --version` !~ m/v1.1.alpha09/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp("Warning : a different version of pplacer was found. Amphora-2 was tested with pplacer v1.1.alpha09\n");
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

	$rapSearch = get_program_path("rapsearch",$Amphora2::Settings::a2_path);
	if($rapSearch eq ""){
	    carp("rapsearch was not found\n");
	    return 1;
	}

	$preRapSearch = get_program_path("prerapsearch",$Amphora2::Settings::a2_path);

	$blastall = get_program_path("blastall",$Amphora2::Settings::a2_path);
	if($blastall eq ""){
	    carp("blastall was not found\n");
	    return 1;
	}
	$formatdb = get_program_path("formatdb",$Amphora2::Settings::a2_path);

	$raxml = get_program_path("raxmlHPC",$Amphora2::Settings::a2_path);
	if($raxml eq ""){
	    carp("raxmlHPC was not found\n");
	    return 1;
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
	# but then again, so is just about every other line of code in this program...
	my $ff = File::Fetch->new(uri=>$url);
	$ff->fetch(to=>"$destination/..");
	debug "URL : $url\n";
	$url =~ /\/(\w+)\.tgz/;
	my $archive = $1;
	debug "ARCHIVE : $archive\n";
	if(-e "$destination/.. "){
	    `rm -rf $destination/..`;
	}
	`cd $destination/../ ; tar xzf $archive.tgz ; touch $archive`;
	`rm $destination/../$archive.tgz`;
}


my $marker_update_url = "http://edhar.genomecenter.ucdavis.edu/~koadman/markers.tgz";
my $ncbi_url = "http://edhar.genomecenter.ucdavis.edu/~koadman/ncbi.tgz";

sub dataChecks {
	$marker_dir = get_data_path( "markers", $Amphora2::Settings::marker_path );
	my ($content_type, $document_length, $modified_time, $expires, $server)= head("$marker_update_url");
	debug "MARKER_PATH : ".$marker_dir."\n";
	if(-x $marker_dir){
	    my $mtime = (stat($marker_dir))[9];
	    debug "TEST LOCAL :".localtime($mtime)."\n";	
	    if(!defined($modified_time)){
		warn "Warning: unable to connect to marker update server, please check your internet connection\n";
	    }elsif($modified_time > $mtime){
		debug "TEST REMOTE:".localtime($modified_time)."\n";
		warn "Found newer version of the marker data\n";
		warn "Downloading from $marker_update_url\n";
		download_data( $marker_update_url, $marker_dir );
	    }
	}else{
	    if(!defined($modified_time)){
		croak "Marker data not found and unable to connect to marker update server, please check your amphora2 configuration and internet connection!\n";
	    }
	    warn "Unable to find marker data!\n";
	    warn "Downloading from $marker_update_url\n";
	    download_data($marker_update_url, $marker_dir);
	}
	$ncbi_dir = get_data_path( "ncbi", $Amphora2::Settings::ncbi_path );
	($content_type, $document_length, $modified_time, $expires, $server)= head("$ncbi_url");
	if( -x $ncbi_dir ){
	    my $ncbi_time =(stat($ncbi_dir))[9];
	    if(!defined($modified_time)){
		warn "Warning: unable to connect to NCBI taxonomy update server, please check your internet connection\n";
	    }elsif($modified_time > $ncbi_time){
		warn "Found newer version of NCBI taxonomy data!\n";
		warn "Downloading from $ncbi_url\n";
		download_data( $ncbi_url, $ncbi_dir );
	    }
	}else{
	    if(!defined($modified_time)){
		croak "NCBI taxonomy data not found and unable to connect to update server, please check your amphora2 configuration and internet connection!\n";
	    }
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
	    my @seq;
	    my @name;
	    my $name;
	    my $curseq = "";
	    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
	    while (<FASTA>) {
		if (/^\s*>\s*(\S+)/) {
		    if(length($curseq)>0){
			push @seq, $curseq;
			push @name, $name;
			$curseq = "";
		    }
		    $name = $1;
		} else {
		    if (/\S/ && !defined $name) {
			warn "Ignoring: $_";
		    } else {
			s/\s//g;
			$curseq .= $_;
		    }
		}
	    }
	    if(length($curseq)>0){
		push @seq, $curseq;
		push @name, $name;
		$curseq = "";
	    }
	    close FASTA;

	# check all seqs are same length
	    my $length;
	    my $lname;
	    for( my $sI=0; $sI<@name; $sI++){
	    	my $nname= $name[$sI];
		my $sseq = $seq[$sI];
		my $l = length $sseq;
		if (defined $length) {
		    croak "Sequences not all same length ($lname is $length, $nname is $l)" unless $length == $l;
		} else {
		    $length = length $sseq;
		    $lname = $nname;
		}
	    }

	# print Stockholm output
	    print STOCKOUT "# STOCKHOLM 1.0\n";
	    for( my $sI=0; $sI<@name; $sI++){
		print STOCKOUT $name[$sI], " ", $seq[$sI], "\n";
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

sub makeDummyFile {
	my $file = shift;
	my $todelete = shift;
	my $gapmultiplier = shift;
	my $stock = $file;	
	$stock =~ s/aln_hmmer3\.trim/trimfinal/g;
	$stock =~ s/\.ffn//g;
	$stock = basename($stock);
	$stock = $Amphora2::Utilities::marker_dir."/$stock";
	open( STOCK, $stock );
	my $burn1 = <STOCK>;
	$burn1 = <STOCK>;
	chomp $burn1;
	my $len = length($burn1);
	my $glen;
	$glen = "P" x $len x $gapmultiplier if $gapmultiplier == 1;
	$glen = "A" x $len x $gapmultiplier if $gapmultiplier == 3;
	my $newseq = Bio::LocatableSeq->new( -seq => $glen, -id => "dummydummydummy", start=>0, end=>($len*$gapmultiplier));
	my $aln = Bio::SimpleAlign->new();
	$aln->add_seq($newseq);
	return $aln;
}

=head2 getAlignmentMarkerFile 

Returns the alignment file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub getAlignemntMarkerFile{
    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
	return "$marker.ali";
    }else{
	return "$marker.updated.fasta";
    }
}
=head2 getFastaMarkerFile

Returns the fasta file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut


sub getFastaMarkerFile{
    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
	return "$marker.faa";
    }else{
	return "$marker.updated.fasta.fasta";
    }
}


=head2 getAlignerOutputFastaAA
Returns the FastA file containing amino acid read or contig alignments to the marker
given by markerName
=cut

sub getAlignerOutputFastaAA{
    my $marker= shift;
    return "$marker.trim.fasta";
}

=head2 getAlignerOutputFastaDNA
Returns the FastA file containing DNA read or contig alignments to the marker
given by markerName
=cut

sub getAlignerOutputFastaDNA{
    my $marker= shift;
    return "$marker.trim.fna.fasta";
}

=head2 getAlignerOutputFastaDNA
Returns the FastA file containing DNA read or contig alignments to the marker
given by markerName
=cut

sub getReadPlacementFile{
    my $marker= shift;
    return "$marker.trim.jplace";
}

=head2 getTrimfinalMarkerFile

Returns the .trimfinal file for the markerName passed in as an argument
If the user chooses the updated markers, the updated filename is returned

=cut

sub getTrimfinalMarkerFile{
    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
        return "$marker.trimfinal";
    }else{
        return "$marker.updated.unique.fasta";
    }
}

=head2 getTrimfinalFastaMarkerFile

Returns the .trimfinal.fasta file for the markerName passed in as an argument
If the use chooses the updated markers, the updated file is returned instead

=cut

sub getTrimfinalFastaMarkerFile{
    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
	return "$marker.trimfinal.fasta";
    }else{
        return $self->{"alignDir"}."/$marker.updated.hmm.fasta";
#        return "$Amphora2::Utilities::marker_dir/$marker.updated.fasta";
    }
}

=head2 getTreeFile

Returns the .final.tre file from the marker directory 
The user chooses the updated or stock version

=cut

sub getTreeMarkerFile{

    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
        return "$marker.final.tre";
    }else{
        return "$marker.updated.tre";
    }

}


=head2 getTreeStatsFile

Return the updated or stock version of the Tree stats file

=cut

sub getTreeStatsMarkerFile{

    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
        return "$marker.in_phyml_stats.txt";
    }else{
        return "$marker.updated.RAxML_info";
    }

}

=head2 getNcbiMapFile

Returns the updated of stock version of the NCBI map file

=cut

sub getNcbiMapFile{

    my $self = shift;
    my $marker= shift;
    if($self->{"updated"} == 0){
        return "$marker.ncbimap";
    }else{
        return "$marker.updated.ncbimap";
    }
}



=head2 concatenateAlignments

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub concatenateAlignments {
	my $outputFasta = shift;
	my $outputMrBayes = shift;
	my $gapmultiplier = shift;	# 1 for protein, 3 for reverse-translated DNA
	my @alignments = @_;
	my $catobj = 0;
	open(MRBAYES, ">$outputMrBayes");
	my $partlist = "partition genes = ".scalar(@alignments).": ";
	my $prevlen = 0;
	my @todelete;
	foreach my $file(@alignments){
		my $aln;
		unless( -e $file ){
			# this marker doesn't exist, need to create a dummy with the right number of gap columns
			$aln = makeDummyFile($file, \@todelete, $gapmultiplier);
		}else{
			my $in  = Bio::AlignIO->new(-file => $file , '-format' => 'fasta');
			unless( $aln = $in->next_aln() ){
				# empty marker alignment file, need to create a dummy with the right number of gap columns
				$aln = makeDummyFile($file, \@todelete, $gapmultiplier);
			}
		}
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
	print MRBAYES "$partlist;\n";
	my $out = Bio::AlignIO->new(-file => ">$outputFasta" , '-format' => 'fasta');
	foreach my $dummyseq( $catobj->each_seq_with_id("dummydummydummy") ){
		$catobj->remove_seq( $dummyseq );
	}
	$out->write_aln($catobj);
	foreach my $delfile(@todelete){
		`rm $delfile`;
	}
}

my %timers;
sub start_timer {
    my $timername = shift;
    my @timerval=localtime(time);
    $timers{$timername} = \@timerval;
    debug join("Before $timername %4d-%02d-%02d %02d:%02d:%02d\n",$timerval[5]+1900,$timerval[4]+1,$timerval[3],$timerval[2],$timerval[1],$timerval[0]);
}

sub end_timer {
    my $timername = shift;
    my @timerval=localtime(time);
    debug join("After $timername %4d-%02d-%02d %02d:%02d:%02d\n",$timerval[5]+1900,$timerval[4]+1,$timerval[3],$timerval[2],$timerval[1],$timerval[0]);
}


=head2 get_sequence_input_type

Checks whether input is either short sequence reads, e.g. < 500nt or assembled fragments
without reading the whole file.

=cut

sub get_sequence_input_type {
	my $file = shift;
	my $maxshortread = 500;
	open(FILE, $file);
	my $filesize = -s "$file";
	my $counter = 0;
	my $allcount = 0;
	my $dnacount = 0;
	my $seqtype="dna";
	my $length="long";
	my $format="unknown";
	
	for( my $i=0; $i<200; $i++ ){
		my $seekpos = int(rand($filesize-100));
		$seekpos = 0 if($i==0); # always start with the first line in case the sequence is on a single line!
		seek(FILE, $seekpos, 0);
		$counter = 0;
		my $line = <FILE>;	# burn a line to be sure we get to sequence
		while( $line = <FILE> ){
			if($line =~ /^>/){
				$format = "fasta";
				last if $i>0;
			}elsif($line =~ /^@/ || $line =~ /^\+/){
				$format = "fastq";
				last if $i>0;
			}else{
				$counter += length($line);		
				$dnacount += $line =~ tr/[ACGTacgt]//;
			}
		}
		$allcount += $counter;
		last if($counter > 500);	# found a long read
	}
	print STDERR "dna frac is ".($dnacount / $allcount)."\n";
	$seqtype = "protein" if ($dnacount < $allcount * 0.75);
	$length = "short" if $counter < 500;
	return ($seqtype, $length, $format);
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
