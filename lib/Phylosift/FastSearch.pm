package Phylosift::FastSearch;
use warnings;
use strict;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::SeqUtils;
use Carp;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(:all);
use File::Basename;
use POSIX qw(ceil floor);
use constant FLANKING_LENGTH => 150;

=head1 NAME

Phylosift::FastSearch - Subroutines to perform fast sequence identity searches between reads and marker genes.
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
my $custom        = "";
my %marker_lookup = ();
my %frames        = ();
my $blastdb_name  = "blastrep.faa";
my $blastp_params = "-p blastp -e 0.1 -b 50000 -v 50000 -m 8";
my $blastn_params = "-p blastn -e 0.1 -b 50000 -v 50000 -m 8";
my %markerLength;

sub RunSearch {
    my %args = @_;
	my $self       = $args{self};
	my $custom     = $args{custom};
	my $markersRef = $args{marker_reference};
	@markers    = @{$markersRef};
	%markerHits = ();
	my $position = rindex( $self->{"readsFile"}, "/" );
	$self->{"readsFile"} =~ m/(\w+)\.?(\w*)$/;
	$readsCore   = $1;
	$isolateMode = $self->{"isolate"};

	# check what kind of input was provided
	my $type = Phylosift::Utilities::get_sequence_input_type( $self->{"readsFile"} );
	$self->{"dna"} = $type->{seqtype} eq "protein" ? 0 : 1;    # Is the input protein sequences?
	debug "Input type is $type->{seqtype}, $type->{format}\n";

	#making sure $type->{paired} is set so we create the appropriate variables
	$type->{paired} = 1 if ( exists $self->{"readsFile_2"} && length( $self->{"readsFile_2"} ) > 0 );

	# need to use BLAST for isolate mode, since RAP only handles very short reads
	# ensure databases and sequences are prepared for search
	debug "before rapPrepandclean\n";
	prep_and_clean( self => $self );
	read_marker_lengths(self=>$self);

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
	my $reads_file      = $args{dir} . "/reads.fasta";
	my $bowtie2_r2_pipe;
	$bowtie2_r2_pipe = $args{dir} . "/bowtie2_r2.pipe" if $args{readtype}->{paired};

	debug "Making fifos\n";
	`mkfifo $blastx_pipe`;
	`mkfifo $bowtie2_r1_pipe`;
	`mkfifo $bowtie2_r2_pipe` if $args{readtype}->{paired};
	`mkfifo $blastp_pipe`;
	`mkfifo $rap_pipe`;
	my @children;

	for ( my $count = 1 ; $count <= 4 ; $count++ ) {
		my $pid = fork();
		if ($pid) {

			# parent process will write sequences below
			push( @children, $pid );
		} elsif ( $pid == 0 ) {
			debug "Launching search process $count\n";

			# child processes will search sequences
			my $hitstream;
			my $candidate_type = ".$count";
			if ( $count == 1 ) {
				$hitstream = lastal_table(self=> $self, query_file=>$blastx_pipe );
				$candidate_type = ".blastx";
			} elsif ( $count == 2 ) {
				$hitstream = bowtie2( self => $self, readtype => $args{readtype}, reads1 => $bowtie2_r1_pipe, reads2 => $bowtie2_r2_pipe );
				$candidate_type = ".rna";
			} elsif ( $count == 3 ) {
				$hitstream = executeBlast( self=>$self, query_file=>$blastp_pipe );
				$candidate_type = ".blastp";
			} elsif ( $count == 4 ) {

				# rapsearch can't read from a pipe
				$hitstream = executeRap( self=>$self, query_file=>$rap_pipe );
				$candidate_type = ".rap";
			}
			my $hitsref;
			if ( $count == 1 ) {
				$hitsref = get_hits_contigs( self => $self, HITSTREAM => $hitstream, searchtype => "lastal" );
			} elsif ( $count == 2 ) {
				$hitsref = get_hits_sam( self => $self, HITSTREAM => $hitstream );
			} elsif ( $count == 4 ) {
				$hitsref = get_hits( self => $self, HITSTREAM => $hitstream, searchtype => "rap" );
				unlink( $self->{"blastDir"} . "/rapjunk.aln" );	# rap leaves some trash lying about
			} elsif ( $args{contigs} ) {
				$hitsref = get_hits_contigs( self => $self, HITSTREAM => $hitstream, searchtype => "lastal" );
			} else {
				$hitsref = get_hits( self => $self, HITSTREAM => $hitstream, searchtype => "blast" );
			}

			# write out sequence regions hitting marker genes to candidate files
			debug "Writing candidates from process $count\n";
			writeCandidates( self => $self, hitsref => $hitsref, searchtype => "$candidate_type", reads => $reads_file );
			exit 0;
		} else {
			croak "couldn't fork: $!\n";
		}
	}
	open( my $RAP_PIPE,        ">$rap_pipe" );
	open( my $BLASTX_PIPE,     ">$blastx_pipe" );
	open( my $BOWTIE2_R1_PIPE, ">$bowtie2_r1_pipe" );
	my $BOWTIE2_R2_PIPE;
	debug "TESTING" . $args{readtype}->{paired};
	open( $BOWTIE2_R2_PIPE, ">$bowtie2_r2_pipe" ) if $args{readtype}->{paired};
	open( my $BLASTP_PIPE,  ">$blastp_pipe" );
	open( my $READS_PIPE,   "+>$reads_file" );

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
					 reads_pipe     => $READS_PIPE,
	);

	# join with children when the searches are done
	foreach (@children) {
		my $tmp = waitpid( $_, 0 );
	}

	# clean up
	`rm -f $blastx_pipe $blastp_pipe $rap_pipe $bowtie2_r1_pipe $reads_file`;
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
	my $READS_PIPE     = $args{reads_pipe};
	my $F1IN           = Phylosift::Utilities::open_sequence_file( file => $args{file1} );
	my $F2IN;
	$F2IN = Phylosift::Utilities::open_sequence_file( file => $args{file2} ) if length( $args{file2} ) > 0;
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
			# send the reads to the reads file in fasta format to write candidates later
			print $READS_PIPE $lines1[0] . $lines1[1];
			print $READS_PIPE $lines2[0] . $lines2[1] if defined($F2IN);

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
			if ( length( $lines1[1] ) > 2000 || ( defined($F2IN) && length( $lines2[1] ) > 2000 ) ) {
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

			#
			# send the reads to the reads file to write candidates later
			print $READS_PIPE $lines1[0] . $lines1[1];
			print $READS_PIPE $lines2[0] . $lines2[1] if defined($F2IN);
			@lines1 = ( $newline1, "" );
			@lines2 = ( $newline2, "" );
		}
	}
	close($RAPSEARCH_PIPE);
	close($BOWTIE2_PIPE1);
	close($BOWTIE2_PIPE2) if defined($F2IN);
	close($BLASTX_PIPE);
	close($BLASTP_PIPE);
	close($RAPSEARCH_PIPE);
	close($READS_PIPE);
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

sub read_marker_lengths {
    my %args = @_;
	my $self = $args{self};
	foreach my $marker (@markers) {
		$markerLength{$marker} = Phylosift::Utilities::get_marker_length( self=>$self, marker=>$marker );
	}
}

=head2 lastal_table

runs lastal with out-of-frame (OOF) detection on a query file (or named pipe)
returns a stream file handle

=cut

sub lastal_table {
    my %args = @_;
	my $self       = $args{self};
	my $query_file = $args{query_file};
	my $lastal_cmd = "$Phylosift::Utilities::lastal -F15 -e300 -f0 $Phylosift::Utilities::marker_dir/replast $query_file |";
	debug "Running $lastal_cmd";
	open( my $hitstream, $lastal_cmd );
	return $hitstream;
}

=head2 blastXoof_table

runs blastx with out-of-frame (OOF) detection on a query file (or named pipe)
returns a stream file handle

=cut

sub blastXoof_table {
    my %args = @_;
	my $self       = $args{self};
	my $query_file = $args{query_file};
	debug "INSIDE tabular OOF blastx\n";
	my $blastxoof_cmd =
	    "$Phylosift::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50000 -v 50000 -d "
	  . Phylosift::Utilities::get_blastp_db( self => $self )
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
	my %args       = @_;
	my $self       = $args{self};
	my $query_file = $args{query};
	debug "INSIDE full OOF blastx\n";
	my $blastxoof_cmd =
	    "$Phylosift::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50 -v 50 -d "
	  . Phylosift::Utilities::get_blastp_db( self => $self ) . " -a "
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
	my %args     = @_;
	my $self     = $args{self};
	my $readtype = $args{readtype};
	debug "INSIDE bowtie2\n";
	my $bowtie2_cmd =
	    "$Phylosift::Utilities::bowtie2align -x "
	  . Phylosift::Utilities::get_bowtie2_db( self => $self )
	  . " --quiet --sam-nohead --sam-nosq --maxins 1000 --mm --local ";
	$bowtie2_cmd .= " -f " if $args{readtype}->{format} eq "fasta";
	if ( $args{readtype}->{paired} ) {
		$bowtie2_cmd .= " -1 $args{reads1} -2 $args{reads2} ";
	} else {
		$bowtie2_cmd .= " -U $args{reads1} ";
	}
	$bowtie2_cmd .= " |";
	debug "Running $bowtie2_cmd";
	open( my $hitstream, $bowtie2_cmd );
	return $hitstream;
}

=head2 translate_frame

=cut

sub translate_frame {
	my %args              = @_;
	my $id                = $args{id};
	my $seq               = $args{seq};
	my $start             = $args{start};
	my $end               = $args{end};
	my $frame             = $args{frame};
	my $marker            = $args{marker};
	my $reverse_translate = $args{reverse_translate};
	my $return_seq        = "";
	my $local_seq         = substr( $seq, $start - 1, $end - $start + 1 );
	my $new_seq           = Bio::LocatableSeq->new( -seq => $local_seq, -id => 'temp' );
	$new_seq = $new_seq->revcom() if ( $frame < 0 );

	if ($reverse_translate) {
		$id = Phylosift::Summarize::treeName($id);
		if ( exists $markerNuc{$marker} ) {
			$markerNuc{$marker} .= ">" . $id . "\n" . $new_seq->seq . "\n";
		} else {
			$markerNuc{$marker} = ">" . $id . "\n" . $new_seq->seq . "\n";
		}
	}
	$return_seq = $new_seq->translate();
	return $return_seq->seq();
}

=head2 executeRap

Launches rapsearch2, returns a stream

=cut

sub executeRap {
    my %args=@_;
	my $self       = $args{self};
	my $query_file = $args{query_file};
	my $dbDir      = "$Phylosift::Utilities::marker_dir/representatives";
	$dbDir = $self->{"blastDir"} if ( $custom ne "" );
	my $out_file      = $self->{"blastDir"} . "/$readsCore.rapSearch";
	my $rapsearch_cmd = "cd "
	  . $self->{"blastDir"}
	  . "; $Phylosift::Utilities::rapSearch -d "
	  . Phylosift::Utilities::get_rapsearch_db( self => $self )
	  . " -o rapjunk -v 20 -b 20 -e -1 -z "
	  . $self->{"threads"}
	  . " < $query_file | ";
	debug "Running $rapsearch_cmd\n";
	open( my $HITSTREAM, $rapsearch_cmd );
	return $HITSTREAM;
}

=head2 executeBlast

Launches blastp, returns a stream

=cut

sub executeBlast {
    my %args=@_;
	my $self       = $args{self};
	my $query_file = $args{query_file};
	my $db         = Phylosift::Utilities::get_blastp_db();
	debug "INSIDE BLAST\n";
	my $blast_cmd = "$Phylosift::Utilities::blastall $blastp_params -i $query_file -d $db -a " . $self->{"threads"} . " |";
	open( my $BLAST_HITS, $blast_cmd );
	return $BLAST_HITS;
}

=head2 get_hits_contigs

parse the blast file

=cut

sub get_hits_contigs {
	my %args       = @_;
	my $self       = $args{self};
	my $HITSTREAM  = $args{HITSTREAM};
	my $searchtype = $args{searchtype};    # can be blastx or lastal

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
		my ( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore );
		if ( $searchtype eq "blastx" ) {
			( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore ) = split( /\t/, $_ );
		} else {
			my @dat = split( /\t/, $_ );
			$bitScore    = $dat[0];
			$subject     = $dat[1];
			$query       = $dat[6];
			$query_start = $dat[7] + 1;
			$query_end   = $query_start + $dat[8] - 1;
			if ( $dat[9] eq "-" ) {

				# reverse strand match
				$query_start = $dat[10] - $dat[7];
				$query_end   = $query_start - $dat[8] + 1;
			}
		}

		# get the marker name
		my @marker = split( /\_\_/, $subject );    # this is soooo ugly
		my $markerName = $marker[0];

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

					# print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end ] if ( $bitScore > $prevhit[1] );
					last;
				}

				# now check the same for reverse-strand hits
				if (    $prevhit[2] > $prevhit[3]
					 && $query_start > $query_end
					 && $prevhit[3] < $query_start - $max_hit_overlap
					 && $query_end + $max_hit_overlap < $prevhit[2] )
				{

					# print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
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
	my %args       = @_;
	my $self       = $args{self};
	my $HITSTREAM  = $args{HITSTREAM};
	my $searchtype = $args{searchtype};
	my %markerTopScores;
	my %topScore = ();
	my %contig_hits;

	# return empty if there is no data
	return \%contig_hits unless defined( fileno $HITSTREAM );
	while (<$HITSTREAM>) {
		chomp($_);
		next if ( $_ =~ /^#/ );
		last if ( $_ =~ /^>>>/ );    # rapsearch end signal
		my ( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore ) = split( /\t/, $_ );
		if ( !defined($subject) ) {
			debug "Undefined subject $_\n";
		}
		my $markerName = get_marker_name( subject => $subject, search_type => $searchtype );

		#parse once to get the top score for each marker (if isolate is ON, assume best hit comes first)
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
	my %args      = @_;
	my $self      = $args{self};
	my $HITSTREAM = $args{HITSTREAM};
	my %markerTopScores;
	my %topScore = ();
	my %contig_hits;

	# return empty if there is no data
	return unless defined($HITSTREAM);
	return \%contig_hits unless defined( fileno $HITSTREAM );
	debug "GETHITSAM\n";
	while (<$HITSTREAM>) {

		#		debug "$_";
		next if ( $_ =~ /^\@/ );
		my @fields = split( /\t/, $_ );
		next if $fields[2] eq "*";    # no hit
		my $marker_name = get_marker_name( subject => $fields[2], search_type => "sam" );
		my $query       = $fields[0];
		my $score       = $fields[4];
		my $cigar       = $fields[5];
		my $qlen        = length( $fields[9] );

		# subtract off soft masking from query length (unaligned portion)
		my $query_lend = 0;
		$query_lend += $1 if $cigar =~ /^(\d+)S/;
		$qlen -= $1 if $cigar =~ /^(\d+)S/;
		$qlen -= $1 if $cigar =~ /(\d+)S$/;
		my $hit_seq = substr( $fields[9], $query_lend, $qlen );

		#add the suffixes back onto the query names if reads are paired
		if ( $fields[1] > 128 ) {

			#greater than 128 the read is the second mate
			$query .= "/2";
		} elsif ( $fields[1] < 128 && $fields[1] > 68 ) {

			#greater than 64 the read is the first mate
			$query .= "/1";
		}

		# flip our coordinates if we're in reverse complement
		# and go back to the start
		$query_lend = length( $fields[9] ) - $query_lend if ( $fields[1] & 0x10 );
		$query_lend = $query_lend - $qlen if ( $fields[1] & 0x10 );
		next if $qlen < 30;    # don't trust anything shorter than 30nt

		# running on short reads, just do one marker per read
		$topScore{$query} = 0 unless exists $topScore{$query};

		#only keep the top hit
		if ( $topScore{$query} <= $score ) {
			$contig_hits{$query} = [ [ $marker_name, $score, $query_lend, $query_lend + $qlen - 1, $hit_seq ] ];
			$topScore{$query} = $score;
		}
	}
	close($HITSTREAM);
	return \%contig_hits;
}

=head2 get_marker_name

Extracts a marker gene name from a blast or rapsearch subject sequence name

=cut

sub get_marker_name {
	my %args        = @_;
	my $subject     = $args{subject};
	my $search_type = $args{search_type};
	my $marker_name = "";
	if ( $search_type eq "blast" ) {
		my @marker = split( /\_/, $subject );
		$marker_name = $marker[$#marker];
	} else {
		my @marker = split( /\_\_/, $subject );
		$marker_name = $marker[0];
	}

	#	debug "Using marker name $markerName";
	return $marker_name;
}

=head2 writeCandidates

write out results

=cut

sub writeCandidates {
	my %args          = @_;
	my $self          = $args{self};
	my $contigHitsRef = $args{hitsref};
	my $type       = $args{searchtype} || ""; # search type -- candidate filenames will have this name embedded, enables parallel output from different programs
	my $reads_file = $args{reads};
	my %contig_hits = %$contigHitsRef;
	debug "ReadsFile:  $self->{\"readsFile\"}" . "\n";
	my $seqin = Phylosift::Utilities::open_SeqIO_object( file => $reads_file );

	while ( my $seq = $seqin->next_seq ) {

		# skip this one if there are no hits
		next unless ( exists $contig_hits{ $seq->id } );
		for ( my $i = 0 ; $i < @{ $contig_hits{ $seq->id } } ; $i++ ) {
			my $cur_hit_ref = $contig_hits{ $seq->id }->[$i];
			my @cur_hit     = @$cur_hit_ref;
			my $markerHit   = $cur_hit[0];
			my $start       = $cur_hit[2];
			my $end         = $cur_hit[3];
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
			$newSeq = $cur_hit[4] if $type =~ /\.rna/;

			#if we're working from DNA then need to translate to protein
			if ( $self->{"dna"} && $type !~ /\.rna/ ) {

				# compute the frame as modulo 3 of start site, reverse strand if end < start
				my $frame = $cur_hit[2] % 3 + 1;
				$frame *= -1 if ( $cur_hit[2] > $cur_hit[3] );
				my $seqlen = abs( $cur_hit[2] - $cur_hit[3] ) + 1;

				# check length again in AA units
				$min_len = $markerLength{$markerHit} < $seq->length / 3 ? $markerLength{$markerHit} : $seq->length / 3;
				next unless ( ( $seqlen / 3 ) / $min_len >= $align_fraction );
				if ( $seqlen % 3 == 0 ) {
					$newSeq = translate_frame(
											   id                => $seq->id,
											   seq               => $seq->seq,
											   start             => $start,
											   end               => $end,
											   frame             => $frame,
											   marker            => $markerHit,
											   reverse_translate => $self->{"dna"}
					);
					$newSeq =~ s/\*/X/g;    # bioperl uses * for stop codons but we want to give X to hmmer later
				} else {
					warn "Search type : $type, alignment length not multiple of 3!  FIXME: need to pull frameshift from full blastx\n";
					next;
				}
			}
			$markerHits{$markerHit} = "" unless defined( $markerHits{$markerHit} );
			$markerHits{$markerHit} .= ">" . $seq->id . "\n" . $newSeq . "\n";
		}
	}

	#write the read+ref_seqs for each markers in the list
	foreach my $marker ( keys %markerHits ) {

		#writing the hits to the candidate file
		my $candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type );
		open( fileOUT, ">$candidate_file" )
		  or croak " Couldn't open $candidate_file for writing\n";
		print fileOUT $markerHits{$marker};
		close(fileOUT);
		if ( $self->{"dna"} && $type !~ /\.rna/ ) {
			$candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type, dna => 1 );
			open( fileOUT, ">$candidate_file" )
			  or croak " Couldn't open $candidate_file for writing\n";
			print fileOUT $markerNuc{$marker} if defined( $markerNuc{$marker} );
			close(fileOUT);
		}
	}
}

=head2 prep_and_clean

=item *

Checks if the directories needed for the blast run and parsing exist
Removes previous blast runs data if they are still in the directories
Generates the blastable database using the marker representatives

=back

=cut

sub prep_and_clean {
	my %args = @_;
	my $self = $args{self};
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

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::blast


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
1;    # End of Phylosift::blast.pm
