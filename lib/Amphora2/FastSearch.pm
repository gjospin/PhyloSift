package Amphora2::FastSearch;
use warnings;
use strict;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::SeqUtils;
use Carp;
use Amphora2::Amphora2;
use Amphora2::Utilities qw(:all);
use File::Basename;
use POSIX qw(ceil floor);
use constant FLANKING_LENGTH => 150;

=head1 NAME

Amphora2::FastSearch - Subroutines to perform fast sequence identity searches between reads and marker genes.
Currently uses either BLAST or RAPsearch.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Run blast on a list of families for a set of Reads
 
 input : Filename with marker list
         Filename for the reads file


 Output : For each marker, create a fasta file with all the reads and reference sequences.
          Storing the files in a directory called Blast_run

 Option : -clean removes the temporary files created to run the blast
          -threaded = #    Runs blast on multiple processors (I haven't see this use more than 1 processor even when specifying more)

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 RunBlast

=cut

my $clean                 = 0;      #option set up, but not used for later
my $isolateMode           = 0;      # set to 1 if running on an isolate assembly instead of raw reads
my $bestHitsBitScoreRange = 30;     # all hits with a bit score within this amount of the best will be used
my $align_fraction        = 0.3;    # at least this amount of min[length(query),length(marker)] must align to be considered a hit
my $pair                  = 0;      #used if using paired FastQ files
my @markers;
my ( %hitsStart, %hitsEnd, %topscore, %hits, %markerHits, %markerNuc ) = ();
my $readsCore;
my $custom           = "";
my %marker_lookup    = ();
my %frames           = ();
my $reverseTranslate = 0;
my $blastdb_name     = "blastrep.faa";
my $blastp_params    = "-p blastp -e 0.1 -b 50000 -v 50000 -m 8";
my $s16db_name       = "16srep.faa";
my $blastn_params    = "-p blastn -e 0.1 -b 50000 -v 50000 -m 8";
my %markerLength;

sub RunSearch {
	my $self       = shift;
	my $custom     = shift;
	my $markersRef = shift;
	@markers    = @{$markersRef};
	%markerHits = ();
	my $position = rindex( $self->{"readsFile"}, "/" );
	$self->{"readsFile"} =~ m/(\w+)\.?(\w*)$/;
	$readsCore        = $1;
	$isolateMode      = $self->{"isolate"};
	$reverseTranslate = $self->{"reverseTranslate"};

	# check what kind of input was provided
	my $type = Amphora2::Utilities::get_sequence_input_type( $self->{"readsFile"} );
	$self->{"dna"} = $type->{seqtype} eq "protein" ? 0 : 1;      # Is the input protein sequences?
	$reverseTranslate = $type->{seqtype} eq "protein" ? 0 : 1;
	debug "Input type is $type->{seqtype}, $type->{format}\n";

	# need to use BLAST for isolate mode, since RAP only handles very short reads
	# ensure databases and sequences are prepared for search
	debug "before rapPrepandclean\n";
	prepAndClean($self);
	readMarkerLengths($self);

	# search reads/contigs against marker database
	my $searchtype = "blast";
	my $contigs    = 0;
	$contigs = 1 if ( defined( $self->{"coverage"} ) && $type->{seqtype} ne "protein" && ( !defined $self->{"isolate"} || $self->{"isolate"} != 1 ) );

	# launch the searches
	launch_searches( self => $self, readtype => $type, dir => $self->{"blastDir"}, contigs => $contigs );

	# delete files that we don't need anymore
	#	cleanup($self);
	return $self;
}

=head2 launch_searches

creates named pipes to stream input to search programs and launches them
after creating pipes, the program forks into separate processes.
the parent process writes sequence data to files, child processes launch a similarity search on that data

=cut

sub launch_searches {
	my %args            = @_;
	my $self            = $args{self};
	my $rap_pipe        = $args{dir} . "/rap.pipe";
	my $blastx_pipe     = $args{dir} . "/blastx.pipe";
	my $blastp_pipe     = $args{dir} . "/blastp.pipe";
	my $bowtie2_r1_pipe = $args{dir} . "/bowtie2_r1.pipe";
	my $bowtie2_r2_pipe;
	$bowtie2_r2_pipe = $args{dir} . "/bowtie2_r2.pipe" if $args{readtype}->{paired};

	#	`mkfifo $rap_pipe`;	# rap doesn't support fifos
	debug "Making fifos\n";
	`mkfifo $blastx_pipe`;
	`mkfifo $bowtie2_r1_pipe`;
	`mkfifo $bowtie2_r2_pipe` if $args{readtype}->{paired};
	`mkfifo $blastp_pipe`;
	my @children;
	for ( my $count = 1 ; $count <= 3 ; $count++ ) {
		my $pid = fork();
		if ($pid) {

			# parent process will write sequences below
			push( @children, $pid );
			debug "Launching search process $count\n";
		} elsif ( $pid == 0 ) {
			debug "Launching search process $count\n";

			# child processes will search sequences
			my $hitstream;
			my $candidate_type = ".$count";
			if ( $count == 1 ) {
				$hitstream = blastXoof_table( $self, $blastx_pipe );
				$candidate_type = ".blastx";
			} elsif ( $count == 2 ) {
				$hitstream = bowtie2( self => $self, readtype => $args{readtype}, reads1 => $bowtie2_r1_pipe, reads2 => $bowtie2_r2_pipe );
				$candidate_type = ".rna";
			} elsif ( $count == 3 ) {
				$hitstream = executeBlast( $self, $blastp_pipe );
				$candidate_type = ".blastp";
			}
			my $hitsref;
			if ( $count == 2 ) {
				$hitsref = get_hits_sam( $self, $hitstream );
			} elsif ( $args{contigs} ) {
				$hitsref = get_hits_contigs( $self, $hitstream, "blast" );
			} else {
				$hitsref = get_hits( $self, $hitstream, "blast" );
			}

			# write out sequence regions hitting marker genes to candidate files
			debug "Writing candidates from process $count\n";
			writeCandidates( $self, $hitsref, "$candidate_type" );
			exit 0;
		} else {
			croak "couldn't fork: $!\n";
		}
	}
	open( my $RAP_PIPE,        ">$rap_pipe" );
	open( my $BLASTX_PIPE,     ">$blastx_pipe" );
	open( my $BOWTIE2_R1_PIPE, ">$bowtie2_r1_pipe" );
	my $BOWTIE2_R2_PIPE;
	open( $BOWTIE2_R2_PIPE, ">$bowtie2_r2_pipe" ) if $args{readtype}->{paired};
	open( my $BLASTP_PIPE, ">$blastp_pipe" );

	# parent process streams out sequences to fifos
	# child processes run the search on incoming sequences
	debug "Demuxing sequences\n";
	demux_sequences(
					 bowtie2_pipe1  => $BOWTIE2_R1_PIPE,
					 bowtie2_pipe2  => $BOWTIE2_R2_PIPE,
					 rapsearch_pipe => $RAP_PIPE,
					 blastx_pipe    => $BLASTX_PIPE,
					 blastp_pipe    => $BLASTP_PIPE,
					 dna            => $self->{"dna"},
					 file1          => $self->{"readsFile"},
					 file2          => $self->{"readsFile_2"},
	);

	# rapsearch can't read from a pipe
	debug "Running rapsearch\n";
	my $rap_hits = executeRap( $self, $rap_pipe );
	my $hitsref = get_hits( $self, $rap_hits, "rap" );
	debug "Done reading rap results, got " . scalar( keys(%$hitsref) ) . " hits\n";
	writeCandidates( $self, $hitsref, ".rap" );

	# join with children when the searches are done
	foreach (@children) {
		my $tmp = waitpid( $_, 0 );
	}

	# clean up
	`rm -f $blastx_pipe $blastp_pipe $rap_pipe $bowtie2_r1_pipe`;
	`rm -f $bowtie2_r2_pipe` if defined($bowtie2_r2_pipe);
}

=head2 demux_sequences

reads a sequence file and streams it out to named pipes

=cut

sub demux_sequences {
	my %args           = @_;
	my $BOWTIE2_PIPE1  = $args{bowtie2_pipe1};
	my $BOWTIE2_PIPE2  = $args{bowtie2_pipe2};
	my $RAPSEARCH_PIPE = $args{rapsearch_pipe};
	my $BLASTX_PIPE    = $args{blastx_pipe};
	my $BLASTP_PIPE    = $args{blastp_pipe};
	my $F1IN           = Amphora2::Utilities::open_sequence_file( file => $args{file1} );
	my $F2IN;
	$F2IN = Amphora2::Utilities::open_sequence_file( file => $args{file2} ) if length( $args{file2} ) > 0;
	my @lines1;
	my @lines2;
	$lines1[0] = <$F1IN>;
	$lines2[0] = <$F2IN> if defined($F2IN);

	while ( defined( $lines1[0] ) ) {
		if ( $lines1[0] =~ /^@/ ) {
			for ( my $i = 1 ; $i < 4 ; $i++ ) {
				$lines1[$i] = <$F1IN>;
				$lines2[$i] = <$F2IN> if defined($F2IN);
			}

			# send the reads to bowtie
			print $BOWTIE2_PIPE1 @lines1;
			print $BOWTIE2_PIPE2 @lines2 if defined($F2IN);

			#
			# send the reads to RAPsearch2 (convert to fasta)
			$lines1[0] =~ s/^@/>/g;
			$lines2[0] =~ s/^@/>/g if defined($F2IN);
			print $RAPSEARCH_PIPE $lines1[0] . $lines1[1];
			print $RAPSEARCH_PIPE $lines2[0] . $lines2[1] if defined($F2IN);

			#
			# prepare for next loop iter
			$lines1[0] = <$F1IN>;
			$lines2[0] = <$F2IN> if defined($F2IN);
		} elsif ( $lines1[0] =~ /^>/ ) {
			my $newline1;
			while ( $newline1 = <$F1IN> ) {
				last if $newline1 =~ /^>/;
				$lines1[1] .= $newline1;
			}
			my $newline2;
			if ( defined($F2IN) ) {
				while ( $newline2 = <$F2IN> ) {
					last if $newline2 =~ /^>/;
					$lines2[1] .= $newline2;
				}
			}

			# send the reads to bowtie
			print $BOWTIE2_PIPE1 @lines1;
			print $BOWTIE2_PIPE2 @lines2 if defined($F2IN);

			# if either read is long, send both to blast
			if ( length( $lines1[1] ) > 500 || ( defined($F2IN) && length( $lines2[1] ) > 500 ) ) {
				if ( $args{dna} ) {
					print $BLASTX_PIPE @lines1;
					print $BLASTX_PIPE @lines2 if defined($F2IN);
				} else {
					print $BLASTP_PIPE @lines1;
					print $BLASTP_PIPE @lines2 if defined($F2IN);
				}
			} else {

				# otherwise send to rap
				print $RAPSEARCH_PIPE $lines1[0] . $lines1[1];
				print $RAPSEARCH_PIPE $lines2[0] . $lines2[1] if defined($F2IN);
			}
			@lines1 = ( $newline1, "" );
			@lines2 = ( $newline2, "" );
		}
	}
	close($RAPSEARCH_PIPE);
	close($BOWTIE2_PIPE1);
	close($BOWTIE2_PIPE2) if $args{readtype}->{paired};
	close($BLASTX_PIPE);
	close($BLASTP_PIPE);
	close($RAPSEARCH_PIPE);
}

sub cleanup {
	my $self = shift;
	`rm -f $self->{"blastDir"}/$readsCore.tabblastx`;
	`rm -f $self->{"blastDir"}/$readsCore.blastx`;
	`rm -f $self->{"blastDir"}/$readsCore.rapsearch.m8`;
	`rm -f $self->{"blastDir"}/$readsCore.rapsearch.aln`;
	`rm -f $self->{"blastDir"}/rep.faa`;
	`rm -f $self->{"blastDir"}/$blastdb_name`;
}

sub readMarkerLengths {
	my $self = shift;
	foreach my $marker (@markers) {
		$markerLength{$marker} = Amphora2::Utilities::get_marker_length( $self, $marker );
	}
}

=head2 blastXoof_table

runs blastx with out-of-frame (OOF) detection on a query file (or named pipe)
returns a stream file handle

=cut

sub blastXoof_table {
	my $self       = shift;
	my $query_file = shift;
	debug "INSIDE tabular OOF blastx\n";
	my $blastxoof_cmd =
	    "$Amphora2::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50000 -v 50000 -d "
	  . Amphora2::Utilities::get_blastp_db()
	  . " -m 8 -a "
	  . $self->{"threads"}
	  . " 2> /dev/null |";
	debug "Running $blastxoof_cmd";
	open( my $hitstream, $blastxoof_cmd );
	return $hitstream;
}

=head2 blastXoof_full

=cut

sub blastXoof_full {
	my $self       = shift;
	my $query_file = shift;
	debug "INSIDE full OOF blastx\n";
	my $blastxoof_cmd =
	    "$Amphora2::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50 -v 50 -d "
	  . Amphora2::Utilities::get_blastp_db() . " -a "
	  . $self->{"threads"}
	  . " 2> /dev/null |";
	debug "Running $blastxoof_cmd";
	open( my $hitstream, $blastxoof_cmd );
	return $hitstream;
}

=head2 bowtie2

runs bowtie2 on a single or pair of query files (or named pipes)
returns a stream file handle

=cut

sub bowtie2 {
	my %args = @_;
	debug "INSIDE bowtie2\n";
	my $bowtie2_cmd =
	  "$Amphora2::Utilities::bowtie2align -x " . Amphora2::Utilities::get_bowtie2_db() . " --quiet --sam-nohead --sam-nosq --maxins 1000 --local ";
	$bowtie2_cmd .= " -f " if $args{readtype}->{format} eq "fasta";
	if ( $args{readtype}->{paired} ) {
		$bowtie2_cmd .= " -1 $args{reads1} -2 $args{reads2} ";
	} else {
		$bowtie2_cmd .= " -U $args{reads1} ";
	}
	$bowtie2_cmd .= "--mm --threads " . $args{self}->{"threads"} . "  |";
	debug "Running $bowtie2_cmd";
	open( my $hitstream, $bowtie2_cmd );
	return $hitstream;
}

=head2 translateFrame

=cut

sub translateFrame {
	my $id               = shift;
	my $seq              = shift;
	my $start            = shift;
	my $end              = shift;
	my $frame            = shift;
	my $marker           = shift;
	my $reverseTranslate = shift;
	my $returnSeq        = "";
	my $localseq         = substr( $seq, $start - 1, $end - $start + 1 );
	my $newseq           = Bio::LocatableSeq->new( -seq => $localseq, -id => 'temp' );
	$newseq = $newseq->revcom() if ( $frame < 0 );

	if ($reverseTranslate) {
		$id = Amphora2::Summarize::treeName($id);
		if ( exists $markerNuc{$marker} ) {
			$markerNuc{$marker} .= ">" . $id . "\n" . $newseq->seq . "\n";
		} else {
			$markerNuc{$marker} = ">" . $id . "\n" . $newseq->seq . "\n";
		}
	}
	$returnSeq = $newseq->translate();
	return $returnSeq->seq();
}

=head2 executeRap

Launches rapsearch2, returns a stream

=cut

sub executeRap {
	my $self       = shift;
	my $query_file = shift;
	my $dbDir      = "$Amphora2::Utilities::marker_dir/representatives";
	$dbDir = $self->{"blastDir"} if ( $custom ne "" );
	my $out_file      = $self->{"blastDir"} . "/$readsCore.rapSearch";
	my $rapsearch_cmd = "cd "
	  . $self->{"blastDir"}
	  . "; $Amphora2::Utilities::rapSearch -q $query_file -d $Amphora2::Utilities::marker_dir/rep -o $out_file -v 20 -b 20 -e -1 -z "
	  . $self->{"threads"} . " > "
	  . $self->{"blastDir"}
	  . "/log.txt";
	debug "Running $rapsearch_cmd\n";
	system($rapsearch_cmd);
	open( my $RAP_HITS, "$out_file.m8" ) || croak "Unable to read from $out_file.m8";
	`rm $out_file.aln`;    # this file is a waste of space
	return $RAP_HITS;
}

=head2 executeBlast

Launches blastp, returns a stream

=cut

sub executeBlast {
	my $self       = shift;
	my $query_file = shift;
	my $db         = Amphora2::Utilities::get_blastp_db();
	debug "INSIDE BLAST\n";
	my $blast_cmd = "$Amphora2::Utilities::blastall $blastp_params -i $query_file -d $db -a " . $self->{"threads"} . " |";
	open( my $BLAST_HITS, $blast_cmd );
	return $BLAST_HITS;
}

=head2 get_hits_contigs

parse the blast file

=cut

sub get_hits_contigs {
	my $self      = shift;
	my $HITSTREAM = shift;

	# key is a contig name
	# value is an array of arrays, each one has [marker,bit_score,left-end,right-end]
	my %contig_hits;
	my %contig_top_bitscore;
	my $max_hit_overlap = 10;

	# return empty if there is no data
	return \%contig_hits unless defined( fileno $HITSTREAM );
	while (<$HITSTREAM>) {

		# read a blast line
		next if ( $_ =~ /^#/ );
		chomp($_);
		my ( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore ) = split( /\t/, $_ );

		# get the marker name
		my @marker = split( /\_/, $subject );
		my $markerName = $marker[$#marker];

		# running on long reads or an assembly
		# allow each region of a sequence to have a top hit
		# do not allow overlap
		if ( defined( $contig_top_bitscore{$query}{$markerName} ) ) {
			my $i = 0;
			for ( ; $i < @{ $contig_hits{$query} } ; $i++ ) {
				my $prevhitref = $contig_hits{$query}->[$i];
				my @prevhit    = @$prevhitref;

				# is there enough overlap to consider these the same?
				# if so, take the new one if it has higher bitscore
				if (    $prevhit[2] < $prevhit[3]
					 && $query_start < $query_end
					 && $prevhit[2] < $query_end - $max_hit_overlap
					 && $query_start + $max_hit_overlap < $prevhit[3] )
				{

					#					print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end ] if ( $bitScore > $prevhit[1] );
					last;
				}

				# now check the same for reverse-strand hits
				if (    $prevhit[2] > $prevhit[3]
					 && $query_start > $query_end
					 && $prevhit[3] < $query_start - $max_hit_overlap
					 && $query_end + $max_hit_overlap < $prevhit[2] )
				{

					#					print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end ] if ( $bitScore > $prevhit[1] );
					last;
				}
			}
			if ( $i == @{ $contig_hits{$query} } ) {

				# no overlap was found, include this hit
				my @hitdata = [ $markerName, $bitScore, $query_start, $query_end ];
				push( @{ $contig_hits{$query} }, @hitdata );
			}
		} elsif ( !defined( $contig_top_bitscore{$query}{$markerName} ) ) {
			my @hitdata = [ $markerName, $bitScore, $query_start, $query_end ];
			push( @{ $contig_hits{$query} }, @hitdata );
			$contig_top_bitscore{$query}{$markerName} = $bitScore;
		}
	}
	return \%contig_hits;
}

=head2 get_hits

parse the blast file, return a hash containing hits to reads

=cut

sub get_hits {
	my $self       = shift;
	my $HITSTREAM  = shift;
	my $searchtype = shift;
	my %markerTopScores;
	my %topScore = ();
	my %contig_hits;

	# return empty if there is no data
	return \%contig_hits unless defined( fileno $HITSTREAM );
	while (<$HITSTREAM>) {
		chomp($_);
		next if ( $_ =~ /^#/ );
		my ( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore ) = split( /\t/, $_ );
		my $markerName = getMarkerName( $subject, $searchtype );

		#parse once to get the top score for each marker (if isolate is ON, parse again to check the bitscore ranges)
		if ( $isolateMode == 1 ) {

			# running on a genome assembly, allow only 1 hit per marker (TOP hit)
			if ( !defined( $markerTopScores{$markerName} ) || $markerTopScores{$markerName} < $bitScore ) {
				$markerTopScores{$markerName} = $bitScore;
			}
			my @hitdata = [ $markerName, $bitScore, $query_start, $query_end ];
			if ( !$self->{"besthit"} && $markerTopScores{$markerName} < $bitScore + $bestHitsBitScoreRange ) {
				push( @{ $contig_hits{$query} }, @hitdata );
			} elsif ( $markerTopScores{$markerName} <= $bitScore ) {
				push( @{ $contig_hits{$query} }, @hitdata );
			}
		} else {

			# running on short reads, just do one marker per read
			$topScore{$query} = 0 unless exists $topScore{$query};

			#only keep the top hit
			if ( $topScore{$query} <= $bitScore ) {
				$contig_hits{$query} = [ [ $markerName, $bitScore, $query_start, $query_end ] ];
				$topScore{$query} = $bitScore;
			}    #else do nothing
		}
	}
	close($HITSTREAM);
	return \%contig_hits;
}

sub get_hits_sam {
	my $self      = shift;
	my $HITSTREAM = shift;
	my %markerTopScores;
	my %topScore = ();
	my %contig_hits;

	# return empty if there is no data
	return unless defined($HITSTREAM);
	return \%contig_hits unless defined( fileno $HITSTREAM );
	while (<$HITSTREAM>) {
		next if ( $_ =~ /^\@/ );
		my @fields = split( /\t/, $_ );
		next if $fields[2] eq "*";    # no hit
		my $markerName = getMarkerName( $fields[2], "sam" );
		my $query      = $fields[0];
		my $score      = $fields[4];
		my $cigar      = $fields[5];
		my $qlen       = length( $fields[9] );

		# subtract off soft masking from query length (unaligned portion)
		my $query_lend = 0;
		$query_lend += $1 if $cigar =~ /^(\d+)S/;
		$qlen -= $1 if $cigar =~ /^(\d+)S/;
		$qlen -= $1 if $cigar =~ /(\d+)S$/;
		my $hit_seq = substr( $fields[9], $query_lend, $qlen );

		# flip our coordinates if we're in reverse complement
		# and go back to the start
		$query_lend = length( $fields[9] ) - $query_lend if ( $fields[1] & 0x10 );
		$query_lend = $query_lend - $qlen if ( $fields[1] & 0x10 );
		next if $qlen < 30;    # don't trust anything shorter than 30nt

		# running on short reads, just do one marker per read
		$topScore{$query} = 0 unless exists $topScore{$query};

		#only keep the top hit
		if ( $topScore{$query} <= $score ) {
			$contig_hits{$query} = [ [ $markerName, $score, $query_lend, $query_lend + $qlen - 1, $hit_seq ] ];
			$topScore{$query} = $score;
		}
	}
	close($HITSTREAM);
	return \%contig_hits;
}

=head2 getMarkerName

Extracts a marker gene name from a blast or rapsearch subject sequence name

=cut

sub getMarkerName {
	my $subject    = shift;
	my $searchtype = shift;
	my $markerName = "";
	if ( $searchtype eq "blast" ) {
		my @marker = split( /\_/, $subject );
		$markerName = $marker[$#marker];
	} else {
		my @marker = split( /\_\_/, $subject );
		$markerName = $marker[0];
	}

	#	debug "Using marker name $markerName";
	return $markerName;
}

=head2 writeCandidates

write out results

=cut

sub writeCandidates {
	my $self          = shift;
	my $contigHitsRef = shift;
	my $type          = shift || "";       # search type -- candidate filenames will have this name embedded, enables parallel output from different programs
	my %contig_hits   = %$contigHitsRef;
	debug "ReadsFile:  $self->{\"readsFile\"}" . "\n";
	my $seqin = new Bio::SeqIO( '-file' => $self->{"readsFile"} );
	while ( my $seq = $seqin->next_seq ) {

		# skip this one if there are no hits
		next unless ( exists $contig_hits{ $seq->id } );
		for ( my $i = 0 ; $i < @{ $contig_hits{ $seq->id } } ; $i++ ) {
			my $curhitref = $contig_hits{ $seq->id }->[$i];
			my @curhit    = @$curhitref;
			my $markerHit = $curhit[0];
			my $start     = $curhit[2];
			my $end       = $curhit[3];
			( $start, $end ) = ( $end, $start ) if ( $start > $end );    # swap if start bigger than end

			# check to ensure hit covers enough of the marker
			# TODO: make this smarter about boundaries, e.g. allow a smaller fraction to hit
			# if it looks like the query seq goes off the marker boundary
			if ( !defined($markerHit) || !defined( $markerLength{$markerHit} ) ) {
				debug "markerHit is $markerHit\n";
				debug $markerLength{$markerHit} . "\n";
			}
			my $min_len = $markerLength{$markerHit} < $seq->length ? $markerLength{$markerHit} : $seq->length;
			next unless ( ( $end - $start ) / $min_len >= $align_fraction );
			$start -= FLANKING_LENGTH;
			$end += FLANKING_LENGTH;

			# ensure flanking region is a multiple of 3 to avoid breaking frame in DNA
			$start = abs($start) % 3 + 1 if ( $start < 0 );
			my $seqLength = length( $seq->seq );
			$end = $end - ceil( ( $end - $seqLength ) / 3 ) * 3 if ( $end >= $seqLength );
			my $newSeq;
			$newSeq = substr( $seq->seq, $start, $end - $start ) unless $type =~ /\.rna/;
			$newSeq = $curhit[4] if $type =~ /\.rna/;

			#if we're working from DNA then need to translate to protein
			if ( $self->{"dna"} && $type !~ /\.rna/ ) {

				# compute the frame as modulo 3 of start site, reverse strand if end < start
				my $frame = $curhit[2] % 3 + 1;
				$frame *= -1 if ( $curhit[2] > $curhit[3] );
				my $seqlen = abs( $curhit[2] - $curhit[3] ) + 1;

				# check length again in AA units
				$min_len = $markerLength{$markerHit} < $seq->length / 3 ? $markerLength{$markerHit} : $seq->length / 3;
				next unless ( ( $seqlen / 3 ) / $min_len >= $align_fraction );
				if ( $seqlen % 3 == 0 ) {
					$newSeq = translateFrame( $seq->id, $seq->seq, $start, $end, $frame, $markerHit, $self->{"dna"} );
					$newSeq =~ s/\*/X/g;    # bioperl uses * for stop codons but we want to give X to hmmer later
				} else {
					warn "Error, alignment length not multiple of 3!  FIXME: need to pull frameshift from full blastx\n";
				}
			}
			$markerHits{$markerHit} = "" unless defined( $markerHits{$markerHit} );
			$markerHits{$markerHit} .= ">" . $seq->id . "\n" . $newSeq . "\n";
		}
	}

	#write the read+ref_seqs for each markers in the list
	foreach my $marker ( keys %markerHits ) {

		#writing the hits to the candidate file
		open( fileOUT, ">" . $self->{"blastDir"} . "/$marker$type.candidate" )
		  or die " Couldn't open " . $self->{"blastDir"} . "/$marker$type.candidate for writing\n";
		print fileOUT $markerHits{$marker};
		close(fileOUT);
		if ( $self->{"dna"} && $type !~ /\.rna/ ) {
			open( fileOUT, ">" . $self->{"blastDir"} . "/$marker$type.candidate.ffn" )
			  or die " Couldn't open " . $self->{"blastDir"} . "/$marker$type.candidate.ffn for writing\n";
			print fileOUT $markerNuc{$marker} if defined( $markerNuc{$marker} );
			close(fileOUT);
		}
	}
}

=head2 prepAndClean

=item *

Checks if the directories needed for the blast run and parsing exist
Removes previous blast runs data if they are still in the directories
Generates the blastable database using the marker representatives

=back

=cut

sub prepAndClean {
	my $self = shift;
	debug "prepclean MARKERS @markers\nTESTING\n ";
	`mkdir $self->{"tempDir"}` unless ( -e $self->{"tempDir"} );

	#create a directory for the Reads file being processed.
	`mkdir $self->{"fileDir"}`  unless ( -e $self->{"fileDir"} );
	`mkdir $self->{"blastDir"}` unless ( -e $self->{"blastDir"} );
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

    perldoc Amphora2::blast


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

1;    # End of Amphora2::blast.pm
