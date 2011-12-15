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
use Amphora2::Utilities qw(:all);


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

my $minAlignedResidues = 20;
my $reverseTranslate;

sub MarkerAlign {

    my $self = shift;
    my $markersRef = shift;
    my @allmarkers = @{$markersRef};
    $reverseTranslate = $self->{"reverseTranslate"};
    debug "beforeDirprepClean @{$markersRef}\n";
    directoryPrepAndClean($self, $markersRef);
    debug "AFTERdirprepclean @{$markersRef}\n";
    my $index =-1;
    
    markerPrepAndRun($self,$markersRef);
    hmmsearchParse($self,$markersRef);
    debug "after HMMSEARCH PARSE\n";
    alignAndMask($self,$reverseTranslate,$markersRef);
    debug "AFTER ALIGN and MASK\n";
#    if($self->{"isolate"} && $self->{"besthit"}){
	    my @markeralignments = getMarkerAlignmentFiles($self,\@allmarkers);
	    Amphora2::Utilities::concatenateAlignments($self->{"alignDir"}."/concat.fasta", $self->{"alignDir"}."/mrbayes.nex", 1, @markeralignments);
	    if($self->{"dna"}){
		for(my $i=0; $i<@markeralignments; $i++){
			$markeralignments[$i] =~ s/trim.fasta/trim.fna.fasta/g;
		}
		Amphora2::Utilities::concatenateAlignments($self->{"alignDir"}."/concat-dna.fasta", $self->{"alignDir"}."/mrbayes-dna.nex", 3, @markeralignments);
	    }
 #   }
    debug "AFTER concatenateALI\n";
    return $self;
}


=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean{
    my $self = shift;
    my $markRef = shift;
    `mkdir -p $self->{"tempDir"}`;
    #create a directory for the Reads file being processed.
    `mkdir -p $self->{"fileDir"}`;
    `mkdir -p $self->{"alignDir"}`;
    for( my $index = 0; $index < @{$markRef}; $index++){
	my $marker = ${$markRef}[$index];
	my $sizer = -s $self->{"blastDir"}."/$marker.candidate"; 
	if(-z $self->{"blastDir"}."/$marker.candidate"){
            warn "WARNING : the candidate file for $marker is empty\n";
            splice @{$markRef}, $index--, 1;
            next;
        }
    }
    return $self;
}

=cut

=head2 markerPrepAndRun

=cut

sub markerPrepAndRun{
    my $self = shift;
    my $markRef = shift;
    debug "ALIGNDIR : ".$self->{"alignDir"}."\n";
    foreach my $marker (@{$markRef}){
	my $trimfinalFile = Amphora2::Utilities::getTrimfinalMarkerFile($self,$marker);
	#converting the marker's reference alignments from Fasta to Stockholm (required by Hmmer3)
	Amphora2::Utilities::fasta2stockholm( "$Amphora2::Utilities::marker_dir/$trimfinalFile", $self->{"alignDir"}."/$marker.seed.stock" );    
	#build the Hmm for the marker using Hmmer3
	if(!-e $self->{"alignDir"}."/$marker.stock.hmm"){
	    `$Amphora2::Utilities::hmmbuild $self->{"alignDir"}/$marker.stock.hmm $self->{"alignDir"}/$marker.seed.stock`;
	}

	if(!-e $self->{"alignDir"}."$marker.hmmsearch.out"){
	    `$Amphora2::Utilities::hmmsearch -E 10 --max --tblout $self->{"alignDir"}/$marker.hmmsearch.tblout $self->{"alignDir"}/$marker.stock.hmm $self->{"blastDir"}/$marker.candidate > $self->{"alignDir"}/$marker.hmmsearch.out`;
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
	open(TBLOUTIN,$self->{"alignDir"}."/$marker.hmmsearch.tblout");
	my $countHits = 0;
	while(<TBLOUTIN>){
	    chomp($_);
	    if($_ =~ m/^(\S+)\s+-\s+(\S+)\s+-\s+(\S+)\s+(\S+)/){
		$countHits++;
		my $hitname = $1;
		my $basehitname = $1;
		my $hitscore = $4;
		if(!defined($hmmScores{$basehitname}) || $hmmScores{$basehitname} < $hitscore ){
		    $hmmScores{$basehitname}=$hitscore;
		    $hmmHits{$basehitname}=$hitname;
		}
	    }
	}
	close(TBLOUTIN);
        
	# added a check if the hmmsearch found hits to prevent the masking and aligning from failing
	if($countHits==0){
	    warn "WARNING : The hmmsearch for $marker found 0 hits, removing marker from the list to process\n";
	    splice @{$markRef}, $index--, 1;
	    next;
	}
	
	open(NEWCANDIDATE,">".$self->{"alignDir"}."/$marker.newCandidate");
	my $seqin = new Bio::SeqIO('-file'=>$self->{"blastDir"}."/$marker.candidate");
	while(my $sequence = $seqin->next_seq){
	    my $baseid = $sequence->id;
	    if(exists $hmmHits{$baseid} && $hmmHits{$baseid} eq $sequence->id){
		print NEWCANDIDATE ">".$sequence->id."\n".$sequence->seq."\n";
	    }
	}
	close(NEWCANDIDATE);
    }
    return $self;
}

=head2 alignAndMask

=cut

sub writeAlignedSeq{
	my $self = shift;
	my $output = shift;
	my $unmaskedout = shift;
	my $prev_name = shift;
	my $prev_seq = shift;
	my $seq_count = shift;
	my $orig_seq = $prev_seq;

	$prev_seq =~ s/[a-z]//g;	# lowercase chars didnt align to model
	$prev_seq =~ s/\.//g;		# shouldnt be any dots
	#skip paralogs if we don't want them
	return if $seq_count > 0 && $self->{"besthit"};
	my $aligned_count = 0;
	$aligned_count++ while $prev_seq =~ m/[A-Z]/g;
	return if $aligned_count < $minAlignedResidues;

	#substitute all the non letter or number characters into _ in the IDs to avoid parsing issues in tree viewing programs or others
	$prev_name = Amphora2::Summarize::treeName($prev_name);
	#add a paralog ID if we're running in isolate mode and more than one good hit
	$prev_name .= "_p$seq_count" if $seq_count > 0 && $self->{"isolate"};
	#print the new trimmed alignment
	print $output ">$prev_name\n$prev_seq\n";
	print $unmaskedout ">$prev_name\n$orig_seq\n" if defined($unmaskedout);
}

use constant CODONSIZE => 3;
my $GAP = '-';
my $CODONGAP = $GAP x CODONSIZE;

=head2 aa_to_dna_aln
Function based on BioPerl's aa_to_dna_aln. This one has been modified to preserve . characters and upper/lower casing of the protein
sequence during reverse translation. Needed to mask out HMM aligned sequences.
=cut

sub aa_to_dna_aln {
    my ($aln,$dnaseqs) = @_;
    unless( defined $aln && 
	    ref($aln) &&
	    $aln->isa('Bio::Align::AlignI') ) { 
	croak('Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature');
    }
    my $alnlen = $aln->length;
    my $dnaalign = Bio::SimpleAlign->new();

    foreach my $seq ( $aln->each_seq ) {    
	my $aa_seqstr = $seq->seq();
	my $id = $seq->display_id;
	my $dnaseq = $dnaseqs->{$id} || $aln->throw("cannot find ".
						     $seq->display_id);
	my $start_offset = ($seq->start - 1) * CODONSIZE;

	$dnaseq = $dnaseq->seq();
	my $dnalen = $dnaseqs->{$id}->length;
	my $nt_seqstr;
	my $j = 0;
	for( my $i = 0; $i < $alnlen; $i++ ) {
	    my $char = substr($aa_seqstr,$i + $start_offset,1);	    
	    if ( $char eq $GAP )  { 
		$nt_seqstr .= $CODONGAP;
	    }elsif ( $char eq "." ){
		$nt_seqstr .= "...";
	    } else {
		if($char eq uc($char)){
			$nt_seqstr .= uc(substr($dnaseq,$j,CODONSIZE));
		}else{
			$nt_seqstr .= lc(substr($dnaseq,$j,CODONSIZE));
		}
		$j += CODONSIZE;
	    }
	}
	$nt_seqstr .= $GAP x (($alnlen * 3) - length($nt_seqstr));

	my $newdna = Bio::LocatableSeq->new(-display_id  => $id,
					   -alphabet    => 'dna',
					   -start       => $start_offset+1,
					   -end         => ($seq->end * 
							    CODONSIZE),
					   -strand      => 1,
					   -seq         => $nt_seqstr);    
	$dnaalign->add_seq($newdna);
    }
    return $dnaalign;
}


sub alignAndMask{
    my $self = shift;
    my $reverseTranslate = shift;
    my $markRef = shift;
    for(my $index=0; $index < @{$markRef}; $index++){
	my $marker=${$markRef}[$index];
	
	my $refcount = 0;
	open( HMM, $self->{"alignDir"}."/$marker.stock.hmm" );
	while( my $line = <HMM> ){
		if($line =~ /NSEQ\s+(\d+)/){
			$refcount = $1;
			last;
		}
	}

	# Align the hits to the reference alignment using Hmmer3
	# pipe in the aligned sequences, trim them further, and write them back out
	my $hmmalign = "$Amphora2::Utilities::hmmalign --outformat afa --mapali ".$self->{"alignDir"}."/$marker.seed.stock ".$self->{"alignDir"}."/$marker.stock.hmm ".$self->{"alignDir"}."/$marker.newCandidate |";
	my $outputFastaAA = $self->{"alignDir"}."/".Amphora2::Utilities::getAlignerOutputFastaAA($marker);
	my $outputFastaDNA = $self->{"alignDir"}."/".Amphora2::Utilities::getAlignerOutputFastaDNA($marker);
	open(my $aliout,">".$outputFastaAA) or die "Couldn't open $outputFastaAA for writing\n";
	open(my $updatedout,">".$self->{"alignDir"}."/$marker.updated.hmm.fasta");
	my $prev_seq;
	my $prev_name;
	my $seqCount = 0;
	open( HMMALIGN, $hmmalign );
	open( my $unmaskedout, ">".$self->{"alignDir"}."/$marker.unmasked");
	my $null;
	while( my $line = <HMMALIGN> ){
		chomp $line;
		if($line =~ /^>(.+)/){
			my $new_name = $1;
			writeAlignedSeq($self, $updatedout, $null, $prev_name, $prev_seq, $seqCount) if $seqCount <= $refcount && $seqCount > 0;
			writeAlignedSeq($self, $aliout, $unmaskedout, $prev_name, $prev_seq, $seqCount) if $seqCount > $refcount;
			$seqCount++;
			$prev_name = $new_name;
			$prev_seq = "";
		}else{
			$prev_seq .= $line;
		}
	}
	writeAlignedSeq($self, $updatedout, $null, $prev_name, $prev_seq, $seqCount) if $seqCount <= $refcount;
	writeAlignedSeq($self, $aliout, $unmaskedout, $prev_name, $prev_seq, $seqCount) if $seqCount > $refcount;
	$seqCount -= $refcount;
	close $unmaskedout;

	# do we need to output a nucleotide alignment in addition to the AA alignment?
	if(-e $self->{"blastDir"}."/$marker.candidate.ffn" && -e $outputFastaAA){
	    #if it exists read the reference nucleotide sequences for the candidates
	    my %referenceNuc = ();
		open(REFSEQSIN,$self->{"blastDir"}."/$marker.candidate.ffn") or die "Couldn't open ".$self->{"alignDir"}."/$marker.candidate.ffn for reading\n";
		my $currID="";
		my $currSeq="";
		while(my $line = <REFSEQSIN>){
		    chomp($line);
		    if($line =~ m/^>(.*)/){
			$currID=$1;
		    }else{
			my $tempseq = Bio::LocatableSeq->new( -seq => $line, -id => $currID);
			$referenceNuc{$currID}=$tempseq;
		    }
		}
		close(REFSEQSIN);

	    open(ALITRANSOUT, ">".$outputFastaDNA) or die "Couldn't open ".$outputFastaDNA/" for writing\n";
	    my $aa_ali = new Bio::AlignIO(-file =>$self->{"alignDir"}."/$marker.unmasked",-format=>'fasta');
	    if(my $aln = $aa_ali->next_aln()){
		    my $dna_ali = &aa_to_dna_aln($aln,\%referenceNuc);
		    foreach my $seq ($dna_ali->each_seq()){
			my $cleanseq = $seq->seq;
			$cleanseq =~ s/\.//g;
			$cleanseq =~ s/[a-z]//g;
			print ALITRANSOUT ">".$seq->id."\n".$cleanseq."\n";
		    }	
	    }
	    close(ALITRANSOUT);
	}

	#checking if sequences were written to the marker alignment file
	if($seqCount ==0){
	    #removing the marker from the list if no sequences were added to the alignment file
	    warn "Masking or hmmsearch thresholds failed, removing $marker from the list\n";
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
	push( @markeralignments, $self->{"alignDir"}."/".Amphora2::Utilities::getAlignerOutputFastaAA($marker) );
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
