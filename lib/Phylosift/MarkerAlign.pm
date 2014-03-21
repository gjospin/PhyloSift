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
use Phylosift::Settings;
use Phylosift::Utilities qw(:all);
use File::Basename;

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::MarkerAlign - Subroutines to align reads to marker HMMs

=head1 VERSION

Version 0.01

=cut

my @search_types = ( "", ".lastal" );
my @search_types_rna = ( "", ".lastal.rna", ".rna" );
set_default_values();

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

sub MarkerAlign {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $markersRef = $args{marker_reference} || miss("marker_reference");
	my $chunk      = $args{chunk};
	if ( defined($chunk) ) {
		my @allmarkers = gather_chunky_markers( self => $self, chunk => $chunk );
		$markersRef = \@allmarkers;
	}
	Phylosift::Utilities::start_step( self => $self, chunk => $chunk, step => "Align" );
	my $completed_chunk = Phylosift::Utilities::has_step_completed( self => $self, chunk => $chunk, step => "Align", force => $Phylosift::Settings::force );
	my $previous_step = Phylosift::Utilities::has_step_completed( self => $self, chunk => $chunk, step => "Search" );
	croak("Previous step for chunk $chunk has did not complete. Aborting\n")
	  unless Phylosift::Utilities::has_step_completed( self => $self, chunk => $chunk, step => "Search" );
	unless ($completed_chunk) {
		directoryPrepAndClean( self => $self, marker_reference => $markersRef, chunk => $chunk );
		my $index = -1;
		markerPrep( self => $self, marker_reference => $markersRef, chunk => $chunk );
		debug "after marker prep\n";
		alignAndMask( self => $self, marker_reference => $markersRef, chunk => $chunk );
		debug "AFTER ALIGN and MASK\n";

		# produce a concatenate alignment for the base marker package
		unless ($Phylosift::Settings::extended) {
			my @markeralignments = get_noncore_alignment_files( self => $self, chunk => $chunk );

			my $outputFastaAA = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "concat", chunk => $chunk );

			#		Phylosift::Utilities::concatenate_alignments(
			concatenate_alignments(
									self           => $self,
									output_fasta   => $outputFastaAA,
									output_bayes   => $self->{"alignDir"}."/mrbayes.nex",
									gap_multiplier => 1,
									alignments     => \@markeralignments
			);

			# now concatenate any DNA alignments
			@markeralignments = get_noncore_alignment_files( self => $self, chunk => $chunk, dna => 1 );
			my $output_fasta_DNA = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "concat", dna => 1, chunk => $chunk );
			concatenate_alignments(
									self           => $self,
									output_fasta   => $output_fasta_DNA,
									output_bayes   => $self->{"alignDir"}."/mrbayes-dna.nex",
									gap_multiplier => 3,
									alignments     => \@markeralignments
			);

			# produce a concatenate with 16s + DNA alignments
			push( @markeralignments, $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "16s_reps_bac", chunk => $chunk ) )
			  if -e $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "16s_reps_bac", chunk => $chunk );
			push( @markeralignments, $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "16s_reps_arc", chunk => $chunk ) )
			  if -e $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "16s_reps_arc", chunk => $chunk );
			push( @markeralignments, $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "18s_reps", chunk => $chunk ) )
			  if -e $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "18s_reps", chunk => $chunk );
			$output_fasta_DNA = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => "concat16", dna => 1, chunk => $chunk );
			concatenate_alignments(
									self           => $self,
									output_fasta   => $output_fasta_DNA,
									output_bayes   => $self->{"alignDir"}."/mrbayes-dna16.nex",
									gap_multiplier => 3,
									alignments     => \@markeralignments
			);
			debug "AFTER concatenateALI\n";
		}
		compute_hits_summary( self => $self, chunk => $chunk );
	}
	Phylosift::Utilities::end_step( self => $self, chunk => $chunk, step => "Align" );
	Phylosift::Utilities::write_step_completion_to_run_info( self => $self, chunk => $chunk, step => "Align" ) unless $completed_chunk;

	# if we're chunking, feed the chunk to the next step
	if ( defined($chunk) && $self->{"mode"} eq "all" ) {
		Phylosift::Utilities::end_timer( name => "runAlign" );
		Phylosift::Utilities::start_timer( name => "runPplacer" );
		Phylosift::pplacer::pplacer( self => $self, marker_reference => $markersRef, chunk => $chunk );
	}
	return $self;
}

=head2 compute_hits_summary

	reads all candidate files and compiles hit numbers for each marker

=cut

sub compute_hits_summary {
	my %args  = @_;
	my $chunk = $args{chunk};
	my $self  = $args{self} || miss("PS Object");
	$chunk = defined $chunk ? ".$chunk" : "";
	my $marker_hits_numbers_ref = Phylosift::Utilities::read_marker_summary( self => $self, path => $self->{"alignDir"} );
	my %marker_hits_numbers     = %{$marker_hits_numbers_ref};
	my @candidates              = glob( $self->{"alignDir"}."/*.updated$chunk.fasta" );
	foreach my $cand_file (@candidates) {
		next if $cand_file =~ m/codon/;
		my $grep_cmd = "grep '>' -c $cand_file";
		$cand_file =~ m/alignDir\/([^\.]+)/;
		my $grep_results = `$grep_cmd`;
		chomp($grep_results);
		$marker_hits_numbers{$1} = 0 unless exists $marker_hits_numbers{$1};
		$marker_hits_numbers{$1} += $grep_results;
	}
	Phylosift::Utilities::print_marker_summary( self => $self, path => $self->{"alignDir"}, summary => \%marker_hits_numbers );
}

=head2 set_default_values

	Set the necessary default values for this module if they haven't been read in by the RC file.

=cut

sub set_default_values {
	my %args = @_;
	my $post = $args{post} || 0;
	if ($post) {
		my $minres = $Phylosift::Settings::isolate && $Phylosift::Settings::besthit ? 40 : 20;
		Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::min_aligned_residues, value => $minres );
	}
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::rna_split_size,       value => 500 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::gap_character,        value => "-" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::cm_align_long_tau,    value => "1e-6" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::cm_align_long_mxsize, value => "2500" );

	#Phylosift::Settings::set_default(parameter=>\$Phylosift::Settings::cm_align_long_ali,value=>"");
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::cm_align_short_tau,    value => "1e-20" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::cm_align_short_mxsize, value => "2500" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::cm_align_short_ali,    value => "-l" );
}

=head2 gather_chunky_markers

=cut

sub gather_chunky_markers {
	my %args              = @_;
	my $self              = $args{self} || miss("PS object");
	my $chunk             = $args{chunk} || miss("Chunk");
	my $type              = $args{type};
	my @candidate_markers = Phylosift::Utilities::get_search_output_all_candidate( self => $self, chunk => $chunk );
	my %unique_markers;
	foreach my $line (@candidate_markers) {
		$line =~ m/\/blastDir\/([^\/\.]+)\.\S+.candidate/;
		my $mark = Phylosift::Utilities::get_marker_fullname( marker => $1 );
		$unique_markers{$mark} = 1 if defined($mark) && length($mark) > 0;
	}
	return keys(%unique_markers);
}

=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $markRef = $args{marker_reference} || miss("marker_reference");
	my $chunk   = $args{chunk};

	#create a directory for the Reads file being processed.
	`mkdir -p "$Phylosift::Settings::file_dir"`;
	`mkdir -p "$self->{"alignDir"}"`;
	for ( my $index = 0; $index < @{$markRef}; $index++ ) {
		my $marker = ${$markRef}[$index];
		my $candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", chunk => $chunk );
		if ( -z $candidate_file ) {
			warn "WARNING : the candidate file for $marker is empty\n";
			splice @{$markRef}, $index--, 1;
			next;
		}
	}
	return $self;
}

sub split_rna_on_size {
	my %args      = @_;
	my $in_fasta  = $args{in_file} || miss("in_file");
	my $short_out = $args{short_out} || miss("short_out");
	my $long_out  = $args{long_out} || miss("long_out");
	my $LONG_OUT;
	my $SHORT_OUT;
	my $seq_in = Phylosift::Utilities::open_SeqIO_object( file => $in_fasta );
	while ( my $seq = $seq_in->next_seq ) {
		my $OUT;
		if ( $seq->length > $Phylosift::Settings::rna_split_size ) {
			$LONG_OUT = ps_open( ">>".$long_out ) unless defined $LONG_OUT && fileno $LONG_OUT;
			$OUT = $LONG_OUT;
		} else {
			$SHORT_OUT = ps_open( ">>".$short_out ) unless defined $SHORT_OUT && fileno $SHORT_OUT;
			$OUT = $SHORT_OUT;
		}
		print $OUT ">".$seq->id."\n".$seq->seq."\n";
	}
}

=cut


=head2 markerPrep

=cut

sub markerPrep {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $markRef = $args{marker_reference} || miss("marker_reference");
	my $chunk   = $args{chunk};

	#debug "Running on ".scalar(@{$markRef})." markers\n";
	foreach my $marker ( @{$markRef} ) {
		unless ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {

			# separate RNA candidates by size
			my $candidate_long  = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => ".rna.long",  chunk => $chunk );
			my $candidate_short = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => ".rna.short", chunk => $chunk );
			unlink($candidate_long);
			unlink($candidate_short);
			foreach my $type (@search_types_rna) {
				my $candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type, chunk => $chunk );
				$candidate = Phylosift::Utilities::escape_char( string => $candidate );
				my @candidate_files = glob("$candidate.*");
				foreach my $cand_file (@candidate_files) {

					#debug "SPLITTING $cand_file\n";
					split_rna_on_size( in_file => $cand_file, short_out => $candidate_short, long_out => $candidate_long );
				}
			}
			next;
		}
		my $new_candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", new => 1, chunk => $chunk );
		unlink($new_candidate);
		foreach my $type (@search_types) {
			my $candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type, chunk => $chunk );
			$candidate = Phylosift::Utilities::escape_char( string => $candidate );
			my @candidate_files = glob("$candidate.*");
			foreach my $cand_file (@candidate_files) {
				next unless -e $cand_file;
				`cat $cand_file >> $new_candidate`;
			}
		}
	}
	return $self;
}

=head2 writeAlignedSeq

=cut

sub writeAlignedSeq {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $OUTPUT      = $args{OUTPUT};
	my $UNMASKEDOUT = $args{UNMASKED_OUT};
	my $prev_name   = $args{prev_name} || miss("prev_name");
	my $prev_seq    = $args{prev_seq} || miss("prev_seq");
	my $seq_count   = $args{seq_count};
	my $orig_seq    = $prev_seq;
	$prev_seq =~ s/[a-z]//g;    # lowercase chars didnt align to model
	$prev_seq =~ s/\.//g;       # shouldnt be any dots
	                            #skip paralogs if we don't want them
	return if $seq_count > 0 && $Phylosift::Settings::besthit;
	my $aligned_count = 0;
	$aligned_count++ while $prev_seq =~ m/[A-Z]/g;
	return if $aligned_count < $Phylosift::Settings::min_aligned_residues;

	#substitute all the non letter or number characters into _ in the IDs to avoid parsing issues in tree viewing programs or others
	my $new_name = $prev_name;
	$self->{"read_names"}{$new_name} = () if ( !exists $self->{"read_names"}{$new_name} );
	push( @{ $self->{"read_names"}{$new_name} }, $prev_name );

	#print the new trimmed alignment
	print $OUTPUT ">$new_name\n$prev_seq\n"      if defined($OUTPUT);
	print $UNMASKEDOUT ">$new_name\n$orig_seq\n" if defined($UNMASKEDOUT);
}
use constant CODONSIZE => 3;

#my $GAP      = $Phylosift::Settings::gap_character;
#my $CODONGAP = $GAP x CODONSIZE;

=head2 aa_to_dna_aln
Function based on BioPerl's aa_to_dna_aln. This one has been modified to preserve . characters and upper/lower casing of the protein
sequence during reverse translation. Needed to mask out HMM aligned sequences.
=cut

sub aa_to_dna_aln {
	my %args = @_;
	my ( $aln, $dnaseqs ) = ( $args{aln}, $args{dna_seqs} );
	my $GAP      = $Phylosift::Settings::gap_character;
	my $CODONGAP = $GAP x CODONSIZE;
	unless (    defined $aln
			 && ref($aln)
			 && $aln->isa('Bio::Align::AlignI') )
	{
		croak(
			'Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature'
		);
	}
	my $alnlen   = $aln->length;
	my $dnaalign = Bio::SimpleAlign->new();
	foreach my $seq ( $aln->each_seq ) {
		my $aa_seqstr = $seq->seq();
		my $id        = $seq->display_id;

		#$id =~ s/\d+\.\d+$//g; # FIXME!! this needs to use lookup table # No longer needed since we use coordinates
		my $dnaseq = $dnaseqs->{$id} || $aln->throw( "cannot find ".$seq->display_id );
		my $start_offset = ( $seq->start - 1 ) * CODONSIZE;
		$dnaseq = $dnaseq->seq();
		my $dnalen = $dnaseqs->{$id}->length;
		my $nt_seqstr;
		my $j = 0;

		for ( my $i = 0; $i < $alnlen; $i++ ) {
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
											 -end           => length($nt_seqstr),
											 -strand        => 1,
											 -seq           => $nt_seqstr,
											 -verbose       => -1,
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
	my $self    = $args{self} || miss("self");
	my $markRef = $args{marker_reference} || miss("marker_reference");
	my $chunk   = $args{chunk};

	my $long_rna = -1;
	for ( my $index = 0; $index < @{$markRef}; $index++ ) {
		my $marker         = ${$markRef}[$index];
		my $refcount       = 0;
		my $stockholm_file = Phylosift::Utilities::get_marker_stockholm_file( self => $self, marker => $marker );
		my @lines;
		my $protein = Phylosift::Utilities::is_protein_marker( marker => $marker );
		$long_rna = ( $long_rna + 1 ) % 2 unless $protein;
		my $candidate;
		if ( !$protein && $long_rna == 0 ) {
			$candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => ".rna.short", chunk => $chunk );
			$index--;
		} elsif ( !$protein && $long_rna == 1 ) {
			$candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => ".rna.long", chunk => $chunk );
		} else {
			$candidate = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => "", new => 1, chunk => $chunk );
		}
		next unless -e $candidate && -s $candidate > 0;

		my $hmm_file = Phylosift::Utilities::get_marker_hmm_file( self => $self, marker => $marker, loc => 1 );
		Phylosift::Utilities::build_hmm( marker => $marker ) unless -e $hmm_file;
		my $HMM = ps_open($hmm_file);
		while ( my $line = <$HMM> ) {
			if ( $line =~ /NSEQ\s+(\d+)/ ) {
				$refcount = $1;
				last;
			}
		}
		if ( $protein || !$long_rna ) {

			# Align the hits to the reference alignment using HMMER 3
			# pipe in the aligned sequences, trim them further, and write them back out
			my $hmmalign = "$Phylosift::Settings::hmmalign --outformat afa --mapali ".$stockholm_file." $hmm_file \"$candidate\" |";
			debug "Running $hmmalign\n";
			my $HMMALIGN = ps_open($hmmalign);
			@lines = <$HMMALIGN>;
		} else {
			debug "Setting up cmalign for marker $marker\n";

			#if the marker is long rna, use infernal instead of hmmalign
			# use tau=1e-6 instead of default 1e-7 to reduce memory consumption to under 4GB
			my $fasta = "";
			$refcount = 0;
			my $cmalign =
			  "$Phylosift::Settings::cmalign -q --dna --mxsize $Phylosift::Settings::cm_align_long_mxsize --tau $Phylosift::Settings::cm_align_long_tau "
			  .Phylosift::Utilities::get_marker_cm_file(
														 self   => $self,
														 marker => $marker
			  )." $candidate | ";
			debug "Running $cmalign\n";
			my $CMALIGN = ps_open($cmalign);
			$fasta .= Phylosift::Utilities::stockholm2fasta( in => $CMALIGN );
			@lines = split( /\n/, $fasta );
			next if @lines == 0;
		}
		my $mbname = Phylosift::Utilities::get_marker_basename( marker => $marker );

		my $short_rna = $protein ? 0 : !$long_rna;
		my $long      = $protein ? 0 : $long_rna;
		my $outputFastaAA =
		  $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => $mbname, chunk => $chunk, long => $long, short => $short_rna );
		my $outputFastaDNA = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => $mbname, dna => 1, chunk => $chunk );
		my $ALIOUT = ps_open( ">".$outputFastaAA );
		my $prev_seq;
		my $prev_name;
		my $seqCount    = 0;
		my $chunky      = defined($chunk) ? ".$chunk" : "";
		my $UNMASKEDOUT = ps_open( ">".$self->{"alignDir"}."/$mbname$chunky.unmasked" );

		foreach my $line (@lines) {
			chomp $line;
			if ( $line =~ /^>(.+)/ ) {
				my $new_name = $1;
				writeAlignedSeq(
								 self      => $self,
								 prev_name => $prev_name,
								 prev_seq  => $prev_seq,
								 seq_count => 0
				  )
				  if $seqCount <= $refcount
					  && $seqCount > 0;

				writeAlignedSeq(
								 self         => $self,
								 OUTPUT       => $ALIOUT,
								 UNMASKED_OUT => $UNMASKEDOUT,
								 prev_name    => $prev_name,
								 prev_seq     => $prev_seq,
								 seq_count    => $seqCount - $refcount - 1
				) if $seqCount > $refcount && $seqCount > 0;
				$seqCount++;
				$prev_name = $new_name;
				$prev_seq  = "";
			} else {
				$prev_seq .= $line;
			}
		}
		writeAlignedSeq( self => $self, prev_name => $prev_name, prev_seq => $prev_seq, seq_count => 0 )
		  if $seqCount <= $refcount && $seqCount > 0;
		writeAlignedSeq(
						 self         => $self,
						 OUTPUT       => $ALIOUT,
						 UNMASKED_OUT => $UNMASKEDOUT,
						 prev_name    => $prev_name,
						 prev_seq     => $prev_seq,
						 seq_count    => $seqCount - $refcount - 1
		) if $seqCount > $refcount;
		$seqCount -= $refcount;
		close $UNMASKEDOUT;
		close $ALIOUT;

		my $type = $self->{readtype};
		unless ( defined($type) ) {
			my $reads_file = Phylosift::Utilities::open_sequence_file( file => $self->{"readsFile"} );
			$type = Phylosift::Utilities::get_sequence_input_type($reads_file) unless defined $type;
		}
		if ( $type->{seqtype} ne "protein" && Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {

			# do we need to output a nucleotide alignment in addition to the AA alignment?
			my %referenceNuc = ();    # this will collect all the nucleotide seqs for the marker by name
			foreach my $type (@search_types) {

				#if it exists read the reference nucleotide sequences for the candidates
				my $core_file_name = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type, dna => 1, chunk => $chunk );
				$core_file_name = Phylosift::Utilities::escape_char( string => $core_file_name );
				my @candidate_files = glob("$core_file_name.*");
				foreach my $cand_file (@candidate_files) {
					if ( -e $cand_file && -e $outputFastaAA ) {
						my $REFSEQSIN = ps_open($cand_file);
						my $currID    = "";
						my $currSeq   = "";
						while ( my $line = <$REFSEQSIN> ) {
							chomp($line);
							if ( $line =~ m/^>(.*)/ ) {
								$currID = $1;
							} else {
								my $tempseq = Bio::PrimarySeq->new( -seq => $line, -id => $currID, -nowarnonempty => 1 );

								#debug "ID : $currID \t BioPerlID : ".$tempseq->id."\n";
								$referenceNuc{$currID} = $tempseq;
							}
						}
						close($REFSEQSIN);
					}
				}
			}
			my $ALITRANSOUT = ps_open( ">>".$outputFastaDNA );
			debug "MARKER : $mbname\n";
			my $aa_ali = new Bio::AlignIO( -file => $self->{"alignDir"}."/$mbname$chunky.unmasked", -format => 'fasta' );
			if ( my $aln = $aa_ali->next_aln() ) {
				my $dna_ali = &aa_to_dna_aln( aln => $aln, dna_seqs => \%referenceNuc );
				foreach my $seq ( $dna_ali->each_seq() ) {
					my $cleanseq = $seq->seq;
					$cleanseq =~ s/\.//g;
					$cleanseq =~ s/[a-z]//g;
					print $ALITRANSOUT ">".$seq->id."\n".$cleanseq."\n";
				}
			}
			close($ALITRANSOUT);
		}

		#checking if sequences were written to the marker alignment file
		if ( $seqCount == 0 ) {

			#removing the marker from the list if no sequences were added to the alignment file
			warn "Masking thresholds failed, removing $marker from the list\n";
			splice @{$markRef}, $index--, 1 if $protein;    # index may have been changed if not protein
		} elsif ( $seqCount > 1 && $Phylosift::Settings::unique ) {
			unlink($outputFastaAA);
			unlink($outputFastaDNA);
			unlink( $self->{"alignDir"}."/$mbname$chunky.unmasked" );
		}

		# check alignments so it merges sequences in case of paired end reads
		if ( defined($Phylosift::Settings::paired) && $Phylosift::Settings::paired ) {

			#debug "PAIRED : " . $Phylosift::Settings::paired . "\n";
			merge_alignment( self => $self, alignment_file => $self->{"alignDir"}."/$mbname$chunky.unmasked", type => 'AA' );
			merge_alignment( self => $self, alignment_file => $outputFastaAA, type => 'AA' );
			if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
				merge_alignment( self => $self, alignment_file => $outputFastaDNA, type => 'DNA' );
			}
		}
	}
}

=head2 merge_alignment

merge alignments by combining sequences from paired end reads.
If aligned columns do not match an X will be used for amino acids and a N will be used for nucleotides
if a residue will always win over a gap
=cut

sub merge_alignment {
	my %args     = @_;
	my $self     = $args{self};
	my $ali_file = $args{alignment_file};
	my $type     = $args{type};
	my %seqs     = ();
	my $seq_IO   = Phylosift::Utilities::open_SeqIO_object( file => $ali_file );
	my %coord    = ();
	while ( my $seq = $seq_IO->next_seq() ) {
		my $core = "";
		if ( $seq->id =~ m/^(\d+)(\.\d+\.\d+\.\d+)$/ ) {
			$core = $1;
			if ( exists $coord{$core} ) {
				$coord{$core} .= $2 if exists $coord{$core};
			} else {
				$coord{$core} = $2;
			}
		}
		if ( exists $seqs{$core} ) {
			my @seq1 = split( //, $seqs{$core} );
			my @seq2 = split( //, $seq->seq );
			my $result_seq = "";
			for ( my $i = 0; $i < length( $seqs{$core} ); $i++ ) {
				if ( $seq1[$i] eq $seq2[$i] ) {
					$result_seq .= $seq1[$i];
				} elsif ( $seq1[$i] =~ /[a-z]/ && $seq2[$i] =~ m/[a-z]/ && $seq1[$i] ne $seq2[$i] ) {
					$result_seq .= $type eq 'AA' ? 'x' : 'n';
				} elsif ( $seq1[$i] =~ /[A-Z]/ && $seq2[$i] =~ m/[A-Z]/ && $seq1[$i] ne $seq2[$i] ) {
					$result_seq .= $type eq 'AA' ? 'X' : 'N';
				} elsif ( $seq1[$i] =~ /[-\.]/ && $seq2[$i] =~ m/[A-Za-z]/ ) {
					$result_seq .= $seq2[$i];
				} elsif ( $seq1[$i] =~ m/[A-Za-z]/ && $seq2[$i] =~ /[-\.]/ ) {
					$result_seq .= $seq1[$i];
				} else {
					debug "FOUND A SPECIAL CASE $seq1[$i] $seq2[$i]\n";
				}
			}
			$seqs{$core} = $result_seq;
		} else {
			$seqs{$core} = $seq->seq;
		}
	}

	#print to the alignment file
	my $FH = ps_open( ">".$ali_file );
	foreach my $core ( keys %seqs ) {
		print $FH ">".$core;
		print $FH $coord{$core} if exists $coord{$core};
		print $FH "\n".$seqs{$core}."\n";
	}
	close($FH);
}

=head get_noncore_alignment_files

Gathers the filenames for alignment files in the AlignDir for a specific chunk and sequence type (DNA/protein)
Returns an Array of filenames
Skips all Core marker names

=cut

sub get_noncore_alignment_files {
	my %args             = @_;
	my $self             = $args{self} || miss("self");
	my $chunk            = $args{chunk};
	my $dna              = $args{dna};
	my @markeralignments = ();
	my @marker_list      = Phylosift::Utilities::gather_markers();
	foreach my $marker (@marker_list) {
		next unless $marker =~ /DNGNGWU/ || $marker =~ /PMPROK/;
		push( @markeralignments, $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta( marker => $marker, chunk => $chunk, dna => $dna ) );
	}
	return @markeralignments;
}

sub concatenate_alignments {
	my %args          = @_;
	my $self          = $args{self};
	my $output_fasta  = $args{output_fasta};
	my $gapmultiplier = $args{gap_multiplier};    # 1 for protein, 3 for reverse-translated DNA
	my $aln_ref       = $args{alignments};
	my %concat_aln;
	my $cur_len   = 0;
	my $seq_count = 0;
	my %id_coords;

	foreach my $alnfile (@$aln_ref) {
		my $marker = basename($alnfile);
		$marker =~ s/\..+//g;                     # FIXME: this should really come from a list of markers
		$gapmultiplier = 1 if ( $marker =~ /16s/ || $marker =~ /18s/ );
		my $len = Phylosift::Utilities::get_marker_length( self => $self, marker => $marker );
		if ( -e $alnfile ) {
			my $ALN = ps_open($alnfile);
			my $id;
			while ( my $line = <$ALN> ) {
				chomp $line;
				if ( $line =~ />(\d+)\.(.+)/ ) {
					$seq_count++;
					$id = $1;
					$id = basename( $self->{"readsFile"} ) if $Phylosift::Settings::isolate && $Phylosift::Settings::besthit;
					my $coords = $2;
					$id_coords{$id} = [] unless defined( $id_coords{$id} );
					push( @{ $id_coords{$id} }, $coords );
				} elsif ( defined($id) ) {
					$concat_aln{$id} = "" unless defined( $concat_aln{$id} );
					my $gapfill = $cur_len - length( $concat_aln{$id} );
					$gapfill = $gapfill < 0 ? 0 : $gapfill;
					$concat_aln{$id} .= "-" x $gapfill;
					$concat_aln{$id} .= $line;
				}
			}
		}
		$cur_len += $len * $gapmultiplier;
	}

	#return if $seq_count <= 0; #don't print an empty file if we don't need to
	# write out the alignment
	my $ALNOUT = ps_open( ">".$output_fasta );
	foreach my $id ( keys(%concat_aln) ) {

		# gapfill for the last marker
		my $gapfill = $cur_len - length( $concat_aln{$id} );
		$gapfill = $gapfill < 0 ? 0 : $gapfill;
		$concat_aln{$id} .= "-" x $gapfill;
		my $gcount = ( $concat_aln{$id} =~ tr/-// );
		next if ( $gcount == length( $concat_aln{$id} ) );    # don't write an all-gap seq. these can slip through sometimes.
		                                                      # write
		print $ALNOUT ">$id.".join( ".", @{ $id_coords{$id} } )."\n$concat_aln{$id}\n";
	}
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
