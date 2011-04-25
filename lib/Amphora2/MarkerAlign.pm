package Amphora2::MarkerAlign;

use Cwd;
use warnings;
use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::SeqIO;
use List::Util qw(min);
use Amphora2::Amphora2;
#use Parallel::ForkManager;


=head1 NAME

Amphora2::MarkerAlign - Subroutines to align reads to marker HMMs

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Run HMMalign for a list of families from .candidate files located in $workingDir/Amph_temp/Blast_run/


input : Filename containing the marker list


output : An alignment file for each marker listed in the input file

Option : -threaded = #    Runs Hmmalign using multiple processors.


Perhaps a little code snippet.

    use Amphora2::Amphora2;

    my $foo = Amphora2::Amphora2->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 MarkerAlign

=cut

my $alnLengthCutoff = 0.4;
my $minAlignedResidues = 20;

sub MarkerAlign {

    my $self = shift;
    my $markersRef = shift;
    print "@{$markersRef}\n";
    directoryPrepAndClean($self, $markersRef);
    print "@{$markersRef}\n";
    my $index =-1;
    
    markerPrepAndRun($self,$markersRef);
    hmmsearchParse($self,$markersRef);
    alignAndMask($self,$markersRef);
    my @markeralignments = getMarkerAlignmentFiles($self,$markersRef);
    Amphora2::Utilities::concatenateAlignments($self->{"alignDir"}."/concat.fasta", $self->{"alignDir"}."/mrbayes.nex", @markeralignments);

    return $self;
}


=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean{
    my $self = shift;
    my $markRef = shift;
    `mkdir $self->{"tempDir"}` unless (-e $self->{"tempDir"});
    #create a directory for the Reads file being processed.
    `mkdir $self->{"fileDir"}` unless (-e $self->{"fileDir"});
    `mkdir $self->{"alignDir"}` unless (-e $self->{"alignDir"});
    for( my $index = 0; $index < @{$markRef}; $index++){
        #    $pm->start and next;                                                                                                                        
	my $marker = ${$markRef}[$index];
	my $sizer = -s $self->{"blastDir"}."/$marker.candidate"; 
	if(-z $self->{"blastDir"}."/$marker.candidate"){
            print STDERR "WARNING : the candidate file for $marker is empty\n";
            splice @{$markRef}, $index--, 1;
            next;
        }
    }
    return $self;
}

=head2 markerPrepAndRun

=cut

sub markerPrepAndRun{
    my $self = shift;
    my $markRef = shift;
    print "ALIGNDIR : ".$self->{"alignDir"}."\n";
    foreach my $marker (@{$markRef}){
	
	#converting the marker's reference alignments from Fasta to Stockholm (required by Hmmer3)
	Amphora2::Utilities::fasta2stockholm( "$Amphora2::Utilities::marker_dir/$marker.trimfinal", $self->{"alignDir"}."/$marker.seed.stock" );    
	#build the Hmm for the marker using Hmmer3
	if(!-e $self->{"alignDir"}."/$marker.stock.hmm"){
	    `$Amphora2::Utilities::hmmbuild $self->{"alignDir"}/$marker.stock.hmm $self->{"alignDir"}/$marker.seed.stock`;
	}
	if(!-e $self->{"alignDir"}."$marker.hmmsearch.out"){
	    `$Amphora2::Utilities::hmmsearch -E 0.1 --tblout $self->{"alignDir"}/$marker.hmmsearch.tblout $self->{"alignDir"}/$marker.stock.hmm $self->{"blastDir"}/$marker.candidate > $self->{"alignDir"}/$marker.hmmsearch.out`;
	}
    }
    return $self;
}

=head2 hmmsearchParse

=cut

sub hmmsearchParse{

    my $self = shift;
    my $markRef = shift;
    for(my $index = 0; $index < @{$markRef}; $index++ ){
	my $marker = ${$markRef}[$index];
	my %hmmHits=();
	my %hmmScores=();
	open(tbloutIN,$self->{"alignDir"}."/$marker.hmmsearch.tblout");
	my $countHits = 0;
	while(<tbloutIN>){
	    chomp($_);
	    if($_ =~ m/^(\S+)\s+-\s+(\S+)\s+-\s+(\S+)\s+(\S+)/){
		$countHits++;
		my $hitname = $1;
		my $basehitname = $1;
		my $hitscore = $4;
		# in case we're using 6-frame translation
		$basehitname =~ s/_[fr][123]$//g;
		if(!defined($hmmScores{$basehitname}) || $hmmScores{$basehitname} < $hitscore ){
		    $hmmScores{$basehitname}=$hitscore;
		    $hmmHits{$basehitname}=$hitname;
		}
	    }
	}
	close(tbloutIN);
        
	# added a check if the hmmsearch found hits to prevent the masking and aligning from failing
	if($countHits==0){
	    print STDERR "WARNING : The hmmsearch for $marker found 0 hits, removing marker from the list to process\n";
	    splice @{$markRef}, $index--, 1;
	    next;
	}
	
	open(newCandidate,">".$self->{"alignDir"}."/$marker.newCandidate");
	my $seqin = new Bio::SeqIO('-file'=>$self->{"blastDir"}."/$marker.candidate");
	while(my $sequence = $seqin->next_seq){
	    my $baseid = $sequence->id;
	    $baseid =~ s/_[fr][123]$//g;
	    if(exists $hmmHits{$baseid} && $hmmHits{$baseid} eq $sequence->id){
		print newCandidate ">".$sequence->id."\n".$sequence->seq."\n";
	    }
	}
	close(newCandidate);
    }
    return $self;
}

=head2 alignAndMask

=cut

sub alignAndMask{
    my $self = shift;
    my $markRef = shift;
    for(my $index=0; $index < @{$markRef}; $index++){
	my $marker=${$markRef}[$index];
	#Align the hits to the reference alignment using Hmmer3
	`$Amphora2::Utilities::hmmalign --trim --outformat afa -o $self->{"alignDir"}/$marker.aln_hmmer3.fasta --mapali $self->{"alignDir"}/$marker.seed.stock $self->{"alignDir"}/$marker.stock.hmm $self->{"alignDir"}/$marker.newCandidate`;
	#find out all the indexes that have a . in the reference sequences
	my $originAli = new Bio::AlignIO(-file=>"$Amphora2::Utilities::marker_dir/$marker.trimfinal", -format=>'fasta');
	my %referenceSeqs = ();
	while(my $aln = $originAli->next_aln()){
	    foreach my $seq($aln->each_seq()){
		$referenceSeqs{$seq->id}=1;
	    }
	}
	my %insertHash=();
	#reading the 
	my $aliIN = new Bio::AlignIO(-file =>$self->{"alignDir"}."/$marker.aln_hmmer3.fasta", -format=>'fasta');
	my $masqseq;
	while(my $aln = $aliIN->next_aln()){
	    # mask out the columns that didn't align to the marker, we'll want to remove
	    # them below
	    foreach my $seq ($aln->each_seq()){
		if(exists $referenceSeqs{$seq->id}){
		    $masqseq = $seq->seq();
		    last;
		}
	    }
	}
	print $masqseq."\n";
	# reading and trimming out non-marker alignment columns from Hmmalign output (Hmmer3)
	my $hmmer3Ali = new Bio::AlignIO(-file =>$self->{"alignDir"}."/$marker.aln_hmmer3.fasta",-format=>'fasta');
	open(aliOUT,">".$self->{"alignDir"}."/$marker.aln_hmmer3.trim")or die "Couldn't open ".$self->{"alignDir"}."/$marker.aln_hmmer3.trim for writting\n";
	my $seqCount = 0;
	while(my $aln = $hmmer3Ali->next_aln()){
	    foreach my $seq ($aln->each_seq()){
		#initialize empty array
		my @characters=();
		#if the sequence isn't a reference sequence remove all the characters at the indices that have a 1 value in the %insertHash
		if(!exists $referenceSeqs{$seq->id}){
		    # strip out all the columns that had a . (indicated by masqseq)
		    my $newSeq = $seq->seq();
		    my $ctr = 0;
		    for(my $i=0; $i<length($masqseq); $i++){
			substr($newSeq, $ctr, 1) = "" if substr($masqseq, $i, 1) eq ".";
			$ctr++ unless substr($masqseq, $i, 1) eq ".";
		    }
		    my $seqLen= 0;
		    #sequence length before masking
		    $seqLen++ while $newSeq =~ m/[^-\.]/g;
                    #change the remaining . into - for pplacer to not complain
		    $newSeq =~ s/\./-/g;
		    my $alignCount=0;
		    $alignCount++ while $newSeq =~ m/[^-]/g;
		    my $collen = length($newSeq);
		    my $minRatio = min (($alignCount / $collen),($alignCount / $seqLen));
                    #print STDERR $seq->id."\t$minRatio\t$alignCount\n";
		    next if ($minRatio < $alnLengthCutoff && $alignCount < $minAlignedResidues);
		    my $newIDs = $seq->id;
		    #get rid of reading frame at this stage
		    $newIDs =~ s/_[fr][012]$//g;
		    #subsitute all the non letter or number characters into _ in the IDs to avoid parsing issues in tree viewing programs or others
		    $newIDs =~ s/[^\w\d]/_/g;
                    #print the new trimmed alignment
		    
		    print aliOUT ">".$newIDs."\n".$newSeq."\n";
		    $seqCount++;
		}
	    }
	}
	close(aliOUT);
	#checking if sequences were written to the marker alignment file
	if($seqCount ==0){
	    #removing the marker from the list if no sequences were added to the alignment file
	    print STDERR "Masking or hmmsearch thresholds failed, removing $marker from the list\n";
	    splice @{$markRef}, $index--, 1;
	}
    }
}


sub getMarkerAlignmentFiles{
    my $self = shift;
    my $markRef = shift;
    my @markeralignments = ();
    for(my $index=0; $index < @{$markRef}; $index++){
	my $marker=${$markRef}[$index];
	push( @markeralignments, $self->{"alignDir"}."/$marker.aln_hmmer3.trim" );
    }
    return @markeralignments;
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

1; # End of Amphora2::MarkerAlign.pm
