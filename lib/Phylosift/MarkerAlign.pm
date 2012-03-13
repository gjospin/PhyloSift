package Phylosift::MarkerAlign;
use Cwd;
use warnings;
use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::SeqIO;
use List::Util qw(min);
use Carp;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(:all);

=head1 NAME

Phylosift::MarkerAlign - Subroutines to align reads to marker HMMs

=head1 VERSION

Version 0.01

=cut
our $VERSION = '0.01';

=head1 SYNOPSIS

Run HMMalign for a list of families from .candidate files located in $workingDir/PS_temp/Blast_run/


input : Filename containing the marker list


output : An alignment file for each marker listed in the input file

Option : -threaded = #    Runs Hmmalign using multiple processors.


Perhaps a little code snippet.

    use Phylosift::Phylosift;

    my $foo = Phylosift::Phylosift->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 MarkerAlign

=cut
my $minAlignedResidues = 20;

sub MarkerAlign {
	my %args       = @_;
	my $self       = $args{self} // miss("self");
	my $markersRef = $args{marker_reference} // miss("marker_reference");
	my @allmarkers = @{$markersRef};
	debug "beforeDirprepClean @{$markersRef}\n";
	directoryPrepAndClean( self => $self, marker_reference => $markersRef );
	debug "AFTERdirprepclean @{$markersRef}\n";
	my $index = -1;
	markerPrepAndRun( self => $self, marker_reference => $markersRef );
	debug "after HMMSEARCH PARSE\n";
	alignAndMask( self => $self, marker_reference => $markersRef );
	debug "AFTER ALIGN and MASK\n";

	# produce a concatenate alignment for the base marker package
	unless ( $self->{"extended"} ) {
		my @markeralignments = getPMPROKMarkerAlignmentFiles( self => $self, marker_reference => \@allmarkers );
		my $outputFastaAA = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => "concat" );
		Phylosift::Utilities::concatenate_alignments(
													  self           => $self,
													  output_fasta   => $outputFastaAA,
													  output_bayes   => $self->{"alignDir"} . "/mrbayes.nex",
													  gap_multiplier => 1,
													  alignments     => \@markeralignments
		);

		# now concatenate any DNA alignments
		for ( my $i = 0 ; $i < @markeralignments ; $i++ ) {
			$markeralignments[$i] =~ s/trim.fasta/trim.fna.fasta/g;
			splice @markeralignments, $i--, 1 unless ( -f $markeralignments[$i] );
		}
		my $outputFastaDNA = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_DNA( marker => "concat" );
		Phylosift::Utilities::concatenate_alignments(
													  self           => $self,
													  output_fasta   => $outputFastaDNA,
													  output_bayes   => $self->{"alignDir"} . "/mrbayes-dna.nex",
													  gap_multiplier => 3,
													  alignments     => \@markeralignments
		);
		debug "AFTER concatenateALI\n";
	}
	return $self;
}

=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean {
	my %args    = @_;
	my $self    = $args{self} // miss("self");
	my $markRef = $args{marker_reference} // miss("marker_reference");
	`mkdir -p $self->{"tempDir"}`;

	#create a directory for the Reads file being processed.
	`mkdir -p $self->{"fileDir"}`;
	`mkdir -p $self->{"alignDir"}`;
	for ( my $index = 0 ; $index < @{$markRef} ; $index++ ) {
		my $marker = ${$markRef}[$index];
		my $candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "" );
		if ( -z $candidate_file ) {
			warn "WARNING : the candidate file for $marker is empty\n";
			splice @{$markRef}, $index--, 1;
			next;
		}
	}
	return $self;
}
my @search_types = ( "", ".blastx", ".blastp", ".rap", ".blast" );

=cut

=head2 markerPrepAndRun

=cut

sub markerPrepAndRun {
	my %args    = @_;
	my $self    = $args{self} // miss("self");
	my $markRef = $args{marker_reference} // miss("marker_reference");
	debug "ALIGNDIR : " . $self->{"alignDir"} . "\n";
	foreach my $marker ( @{$markRef} ) {
		next unless Phylosift::Utilities::is_protein_marker( marker => $marker );
		my $hmm_file = Phylosift::Utilities::get_marker_hmm_file( self => $self, marker => $marker, loc => 1 );
		my $stockholm_file = Phylosift::Utilities::get_marker_stockholm_file( self => $self, marker => $marker );
		unless ( -e $hmm_file && -e $stockholm_file ) {
			my $trimfinalFile = Phylosift::Utilities::get_trimfinal_marker_file( self => $self, marker => $marker );

			#converting the marker's reference alignments from Fasta to Stockholm (required by Hmmer3)
			Phylosift::Utilities::fasta2stockholm( fasta => "$trimfinalFile", output => $stockholm_file );

			#build the Hmm for the marker using Hmmer3
			if ( !-e $hmm_file ) {
				`$Phylosift::Utilities::hmmbuild $hmm_file $stockholm_file`;
			}
		}
		my $new_candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", new => 1 );
		unlink($new_candidate);
		foreach my $type (@search_types) {
			my $candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type );
			next unless -e $candidate;
			my $fifo_out = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_marker_basename( marker => $marker ) . ".tmpout.fifo";
			`mkfifo $fifo_out`;
			system( "$Phylosift::Utilities::hmmsearch -E 10 --cpu " . $self->{"threads"} . " --max --tblout $fifo_out $hmm_file $candidate > /dev/null &" );
			open( my $HMMSEARCH, $fifo_out );
			hmmsearch_parse( self => $self, marker => $marker, type => $type, HMMSEARCH => $HMMSEARCH );
			unlink($fifo_out);
		}
	}
	return $self;
}

=head2 hmmsearchParse

=cut

sub hmmsearch_parse {
	my %args      = @_;
	my $self      = $args{self} // miss("self");
	my $marker    = $args{marker} // miss("marker");
	my $type      = $args{type} // miss("type");
	my $HMMSEARCH = $args{HMMSEARCH} // miss("HMMSEARCH");
	my %hmmHits   = ();
	my %hmmScores = ();
	my $countHits = 0;
	while (<$HMMSEARCH>) {
		chomp($_);
		if ( $_ =~ m/^(\S+)\s+-\s+(\S+)\s+-\s+(\S+)\s+(\S+)/ ) {
			$countHits++;
			my $hitname     = $1;
			my $basehitname = $1;
			my $hitscore    = $4;
			if ( !defined( $hmmScores{$basehitname} ) || $hmmScores{$basehitname} < $hitscore ) {
				$hmmScores{$basehitname} = $hitscore;
				$hmmHits{$basehitname}   = $hitname;
			}
		}
	}
	my $new_candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", new => 1 );
	$new_candidate = ">" . $new_candidate if -f $new_candidate;    # append if the file already exists
	$new_candidate = ">" . $new_candidate;                         # otherwise make a new one
	open( NEWCANDIDATE, $new_candidate ) || croak "Unable to write $new_candidate\n";
	my $candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type );
	my $seqin = Phylosift::Utilities::open_SeqIO_object( file => $candidate );
	while ( my $sequence = $seqin->next_seq ) {
		my $baseid = $sequence->id;
		if ( exists $hmmHits{$baseid} && $hmmHits{$baseid} eq $sequence->id ) {
			print NEWCANDIDATE ">" . $sequence->id . "\n" . $sequence->seq . "\n";
		}
	}
	close(NEWCANDIDATE);
}

=head2 writeAlignedSeq

=cut

sub writeAlignedSeq {
	my %args        = @_;
	my $self        = $args{self} // miss("self");
	my $OUTPUT      = $args{OUTPUT};
	my $UNMASKEDOUT = $args{UNMASKED_OUT};
	my $prev_name   = $args{prev_name} // miss("prev_name");
	my $prev_seq    = $args{prev_seq} // miss("prev_seq");
	my $seq_count   = $args{seq_count} // miss("seq_count");
	my $orig_seq    = $prev_seq;
	if ( !defined($prev_seq) ) {
		print "abc";
	}
	$prev_seq =~ s/[a-z]//g;    # lowercase chars didnt align to model
	$prev_seq =~ s/\.//g;       # shouldnt be any dots
	                            #skip paralogs if we don't want them
	return if $seq_count > 0 && $self->{"besthit"};
	my $aligned_count = 0;
	$aligned_count++ while $prev_seq =~ m/[A-Z]/g;
	return if $aligned_count < $minAlignedResidues;

	#substitute all the non letter or number characters into _ in the IDs to avoid parsing issues in tree viewing programs or others
	$prev_name = Phylosift::Summarize::tree_name( name => $prev_name );

	#add a paralog ID if we're running in isolate mode and more than one good hit
	$prev_name .= "_p$seq_count" if $seq_count > 0 && $self->{"isolate"};

	#print the new trimmed alignment
	print $OUTPUT ">$prev_name\n$prev_seq\n"      if defined($OUTPUT);
	print $UNMASKEDOUT ">$prev_name\n$orig_seq\n" if defined($UNMASKEDOUT);
}
use constant CODONSIZE => 3;
my $GAP      = '-';
my $CODONGAP = $GAP x CODONSIZE;

=head2 aa_to_dna_aln
Function based on BioPerl's aa_to_dna_aln. This one has been modified to preserve . characters and upper/lower casing of the protein
sequence during reverse translation. Needed to mask out HMM aligned sequences.
=cut

sub aa_to_dna_aln {
	my %args = @_;
	my ( $aln, $dnaseqs ) = ( $args{aln}, $args{dna_seqs} );
	unless (    defined $aln
			 && ref($aln)
			 && $aln->isa('Bio::Align::AlignI') )
	{
		croak(
'Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature' );
	}
	my $alnlen   = $aln->length;
	my $dnaalign = Bio::SimpleAlign->new();
	foreach my $seq ( $aln->each_seq ) {
		my $aa_seqstr    = $seq->seq();
		my $id           = $seq->display_id;
		my $dnaseq       = $dnaseqs->{$id} || $aln->throw( "cannot find " . $seq->display_id );
		my $start_offset = ( $seq->start - 1 ) * CODONSIZE;
		$dnaseq = $dnaseq->seq();
		my $dnalen = $dnaseqs->{$id}->length;
		my $nt_seqstr;
		my $j = 0;

		for ( my $i = 0 ; $i < $alnlen ; $i++ ) {
			my $char = substr( $aa_seqstr, $i + $start_offset, 1 );
			if ( $char eq $GAP ) {
				$nt_seqstr .= $CODONGAP;
			} elsif ( $char eq "." ) {
				$nt_seqstr .= "...";
			} else {
				if ( length $dnaseq >= $j + CODONSIZE ) {
					if ( $char eq uc($char) ) {
						$nt_seqstr .= uc( substr( $dnaseq, $j, CODONSIZE ) );
					} else {
						$nt_seqstr .= lc( substr( $dnaseq, $j, CODONSIZE ) );
					}
				}
				$j += CODONSIZE;
			}
		}
		$nt_seqstr .= $GAP x ( ( $alnlen * 3 ) - length($nt_seqstr) );
		my $newdna = Bio::LocatableSeq->new(
											 -display_id    => $id,
											 -alphabet      => 'dna',
											 -start         => 1,
											 -end           => $j,
											 -strand        => 1,
											 -seq           => $nt_seqstr,
											 -nowarnonempty => 1
		);
		$dnaalign->add_seq($newdna);
	}
	return $dnaalign;
}

=head2 alignAndMask

=cut 

sub alignAndMask {
	my %args    = @_;
	my $self    = $args{self} // miss("self");
	my $markRef = $args{marker_reference} // miss("marker_reference");
	for ( my $index = 0 ; $index < @{$markRef} ; $index++ ) {
		my $marker         = ${$markRef}[$index];
		my $refcount       = 0;
		my $stockholm_file = Phylosift::Utilities::get_marker_stockholm_file( self => $self, marker => $marker );
		my $hmmalign       = "";
		my $cmalign        = "";
		if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
			my $new_candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", new => 1 );
			next unless -e $new_candidate && -s $new_candidate > 0;
			my $hmm_file = Phylosift::Utilities::get_marker_hmm_file( self => $self, marker => $marker, loc => 1 );
			open( HMM, $hmm_file );
			while ( my $line = <HMM> ) {
				if ( $line =~ /NSEQ\s+(\d+)/ ) {
					$refcount = $1;
					last;
				}
			}

			# Align the hits to the reference alignment using Hmmer3
			# pipe in the aligned sequences, trim them further, and write them back out
			$hmmalign = "$Phylosift::Utilities::hmmalign --outformat afa --mapali " . $stockholm_file . " $hmm_file $new_candidate |";
		} else {
			debug "Setting up cmalign for marker $marker\n";
			my $candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => ".rna" );
			next unless ( -e $candidate );
			$refcount = Phylosift::Utilities::get_count_from_reps( self => $self, marker => $marker );

			#if the marker is rna, use infernal instead of hmmalign
			# use tau=1e-6 instead of default 1e-7 to reduce memory consumption to under 4GB
			$cmalign =
			    "$Phylosift::Utilities::cmalign -q -l --dna --tau 1e-6 "
			  . Phylosift::Utilities::get_marker_cm_file( self => $self, marker => $marker )
			  . " $candidate | ";
		}
		my $outputFastaAA = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker );
		my $outputFastaDNA = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_DNA( marker => $marker );
		my $mbname = Phylosift::Utilities::get_marker_basename( marker => $marker );
		open( my $aliout, ">" . $outputFastaAA ) or die "Couldn't open $outputFastaAA for writing\n";
		my $updatedout;
		my $prev_seq;
		my $prev_name;
		my $seqCount = 0;
		my @lines;

		if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
			debug "Running $hmmalign\n";
			open( HMMALIGN, $hmmalign );
			@lines = <HMMALIGN>;
		} else {
			debug "Running $cmalign\n";
			open( my $CMALIGN, $cmalign );
			my $sto = Phylosift::Utilities::stockholm2fasta( in => $CMALIGN );
			@lines = split( /\n/, $sto );
			$refcount = 0;
		}
		open( my $UNMASKEDOUT, ">" . $self->{"alignDir"} . "/$mbname.unmasked" );
		my $null;
		foreach my $line (@lines) {
			chomp $line;
			if ( $line =~ /^>(.+)/ ) {
				my $new_name = $1;
				if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
					writeAlignedSeq(
									 self         => $self,
									 OUTPUT       => $updatedout,
									 UNMASKED_OUT => $null,
									 prev_name    => $prev_name,
									 prev_seq     => $prev_seq,
									 seq_count    => 0
					) if $seqCount <= $refcount && $seqCount > 0;
					writeAlignedSeq(
									 self         => $self,
									 OUTPUT       => $aliout,
									 UNMASKED_OUT => $UNMASKEDOUT,
									 prev_name    => $prev_name,
									 prev_seq     => $prev_seq,
									 seq_count    => $seqCount - $refcount - 1
					) if $seqCount > $refcount && $seqCount > 0;
				} else {
					writeAlignedSeq(
									 self         => $self,
									 OUTPUT       => $aliout,
									 UNMASKED_OUT => $UNMASKEDOUT,
									 prev_name    => $prev_name,
									 prev_seq     => $prev_seq,
									 seq_count    => $seqCount
					) if $seqCount > 0;
				}
				$seqCount++;
				$prev_name = $new_name;
				$prev_seq  = "";
			} else {
				$prev_seq .= $line;
			}
		}
		if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
			writeAlignedSeq( self => $self, OUTPUT => $updatedout, UNMASKED_OUT => $null, prev_name => $prev_name, prev_seq => $prev_seq, seq_count => 0 )
			  if $seqCount <= $refcount && $seqCount > 0;
			writeAlignedSeq(
							 self         => $self,
							 OUTPUT       => $aliout,
							 UNMASKED_OUT => $UNMASKEDOUT,
							 prev_name    => $prev_name,
							 prev_seq     => $prev_seq,
							 seq_count    => $seqCount - $refcount - 1
			) if $seqCount > $refcount;
		} else {
			writeAlignedSeq(
							 self         => $self,
							 OUTPUT       => $aliout,
							 UNMASKED_OUT => $UNMASKEDOUT,
							 prev_name    => $prev_name,
							 prev_seq     => $prev_seq,
							 seq_count    => $seqCount
			) if $seqCount > 0;
		}
		$seqCount -= $refcount;
		close $UNMASKEDOUT;

		# do we need to output a nucleotide alignment in addition to the AA alignment?
		foreach my $type (@search_types) {
			if ( -e $self->{"blastDir"} . "/$marker$type.candidate.ffn" && -e $outputFastaAA ) {

				#if it exists read the reference nucleotide sequences for the candidates
				my %referenceNuc = ();
				open( REFSEQSIN, $self->{"blastDir"} . "/$marker$type.candidate.ffn" )
				  or die "Couldn't open " . $self->{"alignDir"} . "/$marker$type.candidate.ffn for reading\n";
				my $currID  = "";
				my $currSeq = "";
				while ( my $line = <REFSEQSIN> ) {
					chomp($line);
					if ( $line =~ m/^>(.*)/ ) {
						$currID = $1;
					} else {
						my $tempseq = Bio::PrimarySeq->new( -seq => $line, -id => $currID, -nowarnonempty => 1 );
						$referenceNuc{$currID} = $tempseq;
					}
				}
				close(REFSEQSIN);
				open( ALITRANSOUT, ">" . $outputFastaDNA ) or die "Couldn't open " . $outputFastaDNA / " for writing\n";
				my $aa_ali = new Bio::AlignIO( -file => $self->{"alignDir"} . "/$marker.unmasked", -format => 'fasta' );
				if ( my $aln = $aa_ali->next_aln() ) {
					my $dna_ali = &aa_to_dna_aln( aln => $aln, dna_seqs => \%referenceNuc );
					foreach my $seq ( $dna_ali->each_seq() ) {
						my $cleanseq = $seq->seq;
						$cleanseq =~ s/\.//g;
						$cleanseq =~ s/[a-z]//g;
						print ALITRANSOUT ">" . $seq->id . "\n" . $cleanseq . "\n";
					}
				}
				close(ALITRANSOUT);
			}
		}

		#checking if sequences were written to the marker alignment file
		if ( $seqCount == 0 ) {

			#removing the marker from the list if no sequences were added to the alignment file
			warn "Masking or hmmsearch thresholds failed, removing $marker from the list\n";
			splice @{$markRef}, $index--, 1;
		}
	}
}

sub getPMPROKMarkerAlignmentFiles {
	my %args             = @_;
	my $self             = $args{self} // miss("self");
	my $markRef          = $args{marker_reference} // miss("marker_reference");
	my @markeralignments = ();
	foreach my $marker ( @{$markRef} ) {
		next unless $marker =~ /PMPROK/;
		push( @markeralignments, $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker ) );
	}
	return @markeralignments;
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
1;    # End of Phylosift::MarkerAlign.pm
