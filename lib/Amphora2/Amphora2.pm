#!/usr/bin/perl
package Amphora2::Amphora2;

use 5.006;
use strict;
use warnings;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use Carp;
use File::Basename;
use Amphora2::Utilities;
use Amphora2::MarkerAlign;
use Amphora2::pplacer;
use Amphora2::Summarize;
use Amphora2::FastSearch;
=head2 new

    Returns : Amphora2 project object
    Args : pair,readsFile(,readsFile_2);

=cut

sub new{
    my $self = {};
    $self->{"fileName"}= undef;
    $self->{"workingDir"}=undef;
    $self->{"mode"} = undef;
    $self->{"readsFile"} = undef;
    $self->{"readsFile_2"} = undef;
    $self->{"tempDir"} = undef;
    $self->{"fileDir"} = undef;
    $self->{"blastDir"} = undef;
    $self->{"alignDir"} = undef;
    $self->{"treeDir"} = undef;
    $self->{"dna"}=undef;
    bless($self);
    return $self;
}

=head2 initialize
    
    Initializes the variables for the Amphora2 object
    Using the standard pathnames and the filename

=cut
    
sub initialize{
    my $self = shift;
    my $mode = shift;
    my $readsFile = shift;
    print "READSFILE\t".$readsFile."\n";
    my $readsFile_2="";
    if(scalar(@_) == 1){
	print "FOUND a second file\n";
	$readsFile_2 = shift;
    }
    
    my $position = rindex($readsFile,"/");
    $self->{"fileName"}= substr($readsFile,$position+1,length($readsFile)-$position-1);
    $self->{"workingDir"}=getcwd;
    $self->{"mode"} = $mode; 
    $self->{"readsFile"} = $readsFile; 
    $self->{"readsFile_2"} = $readsFile_2; 
    $self->{"tempDir"} = $self->{"workingDir"}."/Amph_temp"; 
    $self->{"fileDir"} = $self->{"tempDir"}."/".$self->{"fileName"}; 
    $self->{"blastDir"} = $self->{"fileDir"}."/blastDir"; 
    $self->{"alignDir"} = $self->{"fileDir"}."/alignDir";
    $self->{"treeDir"} = $self->{"fileDir"}."/treeDir"; 
    $self->{"dna"}=0;
    return $self;
    
}

=head2 getReadsFile

    returns the file name for the reads

=cut

sub getReadsFile{
    my $self = shift;
    return $self->{"readsFile"};

}

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
    
    Args : $force - if used at hte same time as the all mode it removes the temp directory if it already exists
           $custom - reads in marker names from a custom list. if "" then use all markers in the marker directory
           $continue - if the mode != all and continue = 1 then finish the pipeline otherwise only execute the step specified
           $isolateMode - allows sequences to hit multiple markers (used when running isolate genomes)

    Runs the Amphora-2 pipeline according to the functions passed as arguments

=cut

my $continue=0;
my ($mode,$readsFile,$readsFile_2, $fileName, $tempDir, $fileDir,$blastDir,$alignDir,$treeDir)="";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=0;
my $workingDir = getcwd;
#where everything will be written when Amphora-2 is running


sub run {
    my $self = shift;
    my $force = shift;
    my $custom = shift;
    my $continue = shift;
    my $reverseTranslate=$self->{"reverseTranslate"};
    print "force : $force\n";
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR  "START : %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    #message to output if the script isn't called properly
    my $usage = qq~
	Usage: $0 <mode> <options> <reads_file>

	~;
    my $usage2 = qq~
        Usage: $0 <mode> <options> -paired <reads_file_1> <reads_file_2>

	~;

    $self->readAmphora2Config();
    $self->runProgCheck();
    Amphora2::Utilities::dataChecks($self);
    $self->fileCheck();
    $self->directoryPrep($force);

    $self->{"readsFile"} = $self->prepIsolateFiles($self->{"readsFile"}) if $self->{"isolate"}==1;

    #create a file with a list of markers called markers.list
    print "CUSTOM = ".$custom."\n";

    my @markers = $self->markerGather($custom);
    print "@markers\n";
    print "MODE :: ".$self->{"mode"}."\n";
    if($self->{"mode"} eq 'blast' || $self->{"mode"} eq 'all'){
	my $searchtype = "rap";
	$searchtype = "blast" if defined($self->{"isolate"}) && $self->{"isolate"} ne "0";
	print "Search type is $searchtype\n";
	# need to use BLAST for isolate mode, since RAP only handles very short reads
	$self=$self->runSearch($continue,$custom,$searchtype,\@markers);
	print "MODE :: ".$self->{"mode"}."\n";
    }

    print "MODE :: ".$self->{"mode"}."\n";
    if($self->{"mode"} eq 'align' || $self->{"mode"} eq 'all'){
	$self=$self->runMarkerAlign($continue,\@markers);
    }
    if($self->{"mode"} eq 'placer' || $self->{"mode"} eq 'all'){
	$self=$self->runPplacer($continue,\@markers);
    }
    if($self->{"mode"} eq 'summary' || $self->{"mode"} eq 'all'){
	$self=$self->taxonomyAssignments();
    }
}

=head2 function2

    Reads the Amphora configuration file and assigns the file paths to the required directories

=cut

# reads the Amphora2 configuration file
sub readAmphora2Config {
    my $self = shift;
    # first get the install prefix of this script
    my $scriptpath = dirname($0);
    # try first a config in the script dir, in case we're running from
    # a dev directory.  then config in system dir, then user's home.
    # let each one override its predecessor.
    { package Amphora2::Settings; do "$scriptpath/amphora2rc"; do "$scriptpath/../etc/amphora2rc"; do "$ENV{HOME}/.amphora2rc" }
    return $self;
}

=head2 fileCheck

    Checks if the files passed to the Amphora object exist and are not empty

=cut

sub fileCheck{
    my $self = shift;
    if(!-e $self->{"readsFile"} ){
	die $self->{"readsFile"}."  was not found \n";
    }
    #check if the input file is a file and not a directory
    if( !-f $self->{"readsFile"}){
	die $self->{"readsFile"}. " is not a plain file, could be a directory\n";
    }
    if($self->{"readsFile_2"} ne ""){
	#check the input file exists
	if(!-e $self->{"readsFile_2"}){
	    die $self->{"readsFile_2"}." was not found\n";
	}
	#check if the input file is a file and not a directory
	if( !-e $self->{"readsFile_2"}){
	    die $self->{"readsFile_2"}." is not a plain file, could be a directory\n";
	}
    }
    return $self;
}

=head2 runProgCheck

    Runs a check on the programs that will be used through the pipeline to make sure they are
    available to the user and are the versions Amphora was tested with.

=cut

sub runProgCheck{
    my $self = shift;
    #check if the various programs used in this pipeline are installed on the machine
    my $progCheck = Amphora2::Utilities::programChecks($self);
    if($progCheck!=0){
	print STDERR "A required program was not found during the checks aborting\n";
	exit();
    }elsif($progCheck==0){
	print STDERR "All systems are good to go, continuing the screening\n";
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
	open( OUTFILE, ">".$self->{"fileDir"}."/isolates.fasta" );
	while(my $file = shift){
		open( ISOLATEFILE, $file ) || croak("Unable to read $file\n");
		my $fname = $self->{"fileDir"}."/".basename($file);
		print "Operating on isolate file $fname\n";
		print OUTFILE ">".basename($file)."\n";
		while( my $line = <ISOLATEFILE> ){
			next if $line =~ /^>/;
			print OUTFILE $line;
		}
		print "\n";
		close ISOLATEFILE;
	}
	close OUTFILE;
	$self->{"readsFile"} = "isolates.fasta";
	return $self->{"fileDir"}."/isolates.fasta";
}

=head2 markerGather

=item *

    ARGS : $markerFile - file to read the marker names from that will be used in the pipeline.
    Reads markers from a file if specified by the user OR gathers all the markers from the markers directory
    The Marker names are stored in an Array
    If the filename is empty, use all the default markers.

=back

=cut

sub markerGather {
    my $self = shift;
    my $markerFile = shift;
    my @marks=();
    #create a file with a list of markers called markers.list
    if($markerFile ne ""){
	#gather a custom list of makers, list convention is 1 marker per line
	open(markersIN,$markerFile);
	while(<markersIN>){
	    chomp($_);
	    push(@marks,$_);
	}
	close(markersIN);
    }else{
	#gather all markers
	my @files = <$Amphora2::Utilities::marker_dir/*.faa>;
	foreach my $file (@files){
	    $file =~ m/\/(\w+).faa/;
	    push(@marks,$1);
	}
    }
    return @marks;
}

=head2 directoryPrep

    Prepares the temporary Amphora directory by deleting old runs and/or creating the correct directory structure
    
=cut

sub directoryPrep {
    my $self = shift; 
    my $force = shift;
#    print "FORCE DIRPREP   $force\t mode   ".$self->{"mode"}."\n";
#    exit;
    #remove the directory from a previous run
    if($force && $self->{"mode"} eq 'all'){
	print STDERR "deleting an old run\n";
	my $dir = $self->{"fileDir"};
	`rm -rf $dir`;
    }elsif(-e $self->{"fileDir"} && $self->{"mode"} eq 'all'){
	print STDERR "A previous run was found using the same file name aborting the current run\n";
	print STDERR "Either delete that run from ".$self->{"fileDir"}.", or force overwrite with the -f command-line option\n";
	exit;
    }
    #check if the temporary directory exists, if it doesn't create it.
    `mkdir $self->{"tempDir"}` unless (-e $self->{"tempDir"});
    #create a directory for the Reads file being processed.
    `mkdir $self->{"fileDir"}` unless (-e $self->{"fileDir"});
    `mkdir $self->{"blastDir"}` unless (-e $self->{"blastDir"});
    `mkdir $self->{"alignDir"}` unless (-e $self->{"alignDir"});
    `mkdir $self->{"treeDir"}` unless (-e $self->{"treeDir"});
    return $self;
}

=head2 taxonomyAssignments

    performs the appropriate checks before running the taxonomy classification parts of the pipeline

=cut

sub taxonomyAssignments {
    my $self = shift;
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "Before taxonomy assignments %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    # Taxonomy assignemnts
    print STDERR "NEIGHOR TEST \t$Amphora2::Utilities::printneighbor\n";
    `$Amphora2::Utilities::Rscript $Amphora2::Utilities::printneighbor $self->{"treeDir"}/*.num.tre > $self->{"fileDir"}/neighbortaxa.txt`;
    Amphora2::Summarize::summarize( $self );
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "After taxonomy assignments %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
}

=head2 runPplacer

    Performs the appropriate checks before runing the Read placement parts of the pipeline

    if -continue or all mode are used, then run the next part of the pipeline

=cut

sub runPplacer{
    my $self = shift;
    my $continue = shift;
    my $markListRef = shift;
    print "PPLACER MARKS @{$markListRef}\n";
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "Starting runPPlacer %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    my $treeDir =$self->{"treeDir"};
    `rm $treeDir/*` if (<$treeDir/*>);
    Amphora2::pplacer::pplacer($self,$markListRef);
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "After runPPlacer %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    if($continue != 0){
	$self->{"mode"} = 'summary';
    }
    return $self;
}

=head2 runMarkerAlign

    Run the Hmm hit verification and the hit alignments to prep the hits for Pplacer
    if -continue or all mode are used, then run the next part of the pipeline

=cut
    
sub runMarkerAlign{
    my $self = shift;
    my $continue = shift;
    my $markRef = shift;
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "Before Alignments for Markers %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    #clearing the alignment directory if needed
    my $alignDir = $self->{"alignDir"};
    `rm $alignDir/*` if(<$alignDir/*>);
    #Align Markers
    my $threadNum=1;
    Amphora2::MarkerAlign::MarkerAlign( $self, $markRef );
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "After Alignments %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    if($continue != 0){
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
    my $self = shift;
    my $continue = shift;
    my $custom = shift;
    my $type = shift;
    my $markerListRef = shift;
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "Before runBlast %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    #clearing the blast directory
    my $blastDir = $self->{"blastDir"};
    `rm $self->{"blastDir"}/*` if(<$blastDir/*>);
    #run Blast
    Amphora2::FastSearch::RunSearch($self,$custom,$type,$markerListRef);
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    printf STDERR "After runBlast %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
    
    if($continue != 0){
	$self->{"mode"} = 'align';
    }
    return $self;
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
