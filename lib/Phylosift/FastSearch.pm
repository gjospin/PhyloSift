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
my $bestHitsBitScoreRange  = 30;     # all hits with a bit score within this amount of the best will be used
my $align_fraction         = 0.5;    # at least this amount of min[length(query),length(marker)] must align to be considered a hit
my $align_fraction_isolate = 0.8;    # use this align_fraction when in isolate mode on long sequences
my $pair                   = 0;      #used if using paired FastQ files
my @markers;
my %markerNuc = ();
my $readsCore;
my $blastdb_name  = "blastrep.faa";
my $blastp_params = "-p blastp -e 0.1 -b 50000 -v 50000 -m 8";
my $blastn_params = "-p blastn -e 0.1 -b 50000 -v 50000 -m 8";
my %markerLength;

sub RunSearch {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $custom     = $args{custom};
	my $markersRef = $args{marker_reference} || miss("marker_reference");
	@markers = @{$markersRef};
	my $position = rindex( $self->{"readsFile"}, "/" );
	$self->{"readsFile"} =~ m/(\w+)\.?(\w*)$/;
	$readsCore = $1;

	# set align_fraction appropriately
	$align_fraction = $align_fraction_isolate if ( $self->{"isolate"} );

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
	read_marker_lengths( self => $self );

	# search reads/contigs against marker database
	my $searchtype = "blast";
	my $contigs    = 0;
	$contigs = 1 if ( defined( $self->{"coverage"} ) && $type->{seqtype} ne "protein" && ( !defined $self->{"isolate"} || $self->{"isolate"} != 1 ) );

	# launch the searches
	launch_searches( self => $self, readtype => $type, dir => $self->{"blastDir"}, contigs => $contigs );

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
	my $dir             = $args{dir} || miss("dir");
	my $bowtie2_r1_pipe = $args{dir} . "/bowtie2_r1.pipe";
	my $reads_file      = $args{dir} . "/reads.fasta";
	my $bowtie2_r2_pipe;
	$bowtie2_r2_pipe = $args{dir} . "/bowtie2_r2.pipe" if $args{readtype}->{paired};
	debug "Making fifos\n";
	
	my @last_pipe_array       = ();
	for(my $i =0 ; $i <= $self->{"threads"}-1;$i++){
		push(@last_pipe_array,$args{dir} . "/last_$i.pipe");
		`mkfifo "$args{dir}/last_$i.pipe"`;
	}
	
	`mkfifo $bowtie2_r1_pipe`;
	`mkfifo $bowtie2_r2_pipe` if $args{readtype}->{paired};
	my @children;

	for ( my $count = 1 ; $count <= $self->{"threads"} + 1; $count++ ) {
		my $pid = fork();
		if ($pid) {

			# parent process will write sequences below
			push( @children, $pid );
		} elsif ( $pid == 0 ) {
			debug "Launching search process $count\n";

			# child processes will search sequences
			my $hitstream;
			my $candidate_type = ".$count";
			if ( $count <= $self->{"threads"} ) {
				$hitstream = lastal_table( self => $self, query_file => $last_pipe_array[$count-1] );
				$candidate_type = ".lastal";
			} elsif ( $count > $self->{"threads"} ) {

				#exit the thread if the bowtie DB does not exist
				exit 0 unless -e Phylosift::Utilities::get_bowtie2_db();
				$hitstream = bowtie2( self => $self, readtype => $args{readtype}, reads1 => $bowtie2_r1_pipe, reads2 => $bowtie2_r2_pipe );
				$candidate_type = ".rna";
			}
			my $hitsref;
			if ( $count <= $self->{"threads"} ) {
				$hitsref = get_hits_contigs( self => $self, HITSTREAM => $hitstream, searchtype => "lastal" );
			} elsif ( $count > $self->{"threads"} ) {
				$hitsref = get_hits_sam( self => $self, HITSTREAM => $hitstream );

				#$hitsref = get_hits_contigs( self => $self, HITSTREAM => $hitstream, searchtype => "lastal" );
			}

			# write out sequence regions hitting marker genes to candidate files
			debug "Writing candidates from process $count\n";
			write_candidates( self => $self, hitsref => $hitsref, searchtype => "$candidate_type", reads => $reads_file, process_id => $count );
			exit 0;
		} else {
			croak "couldn't fork: $!\n";
		}
	}
	my @LAST_PIPE_ARRAY = ();
	foreach my $last_pipe (@last_pipe_array){
		push(@LAST_PIPE_ARRAY, ps_open( ">$last_pipe" ));
	}
	my $BOWTIE2_R1_PIPE = ps_open( ">$bowtie2_r1_pipe" );
	my $BOWTIE2_R2_PIPE;
	debug "TESTING" . $args{readtype}->{paired};
	$BOWTIE2_R2_PIPE = ps_open( ">$bowtie2_r2_pipe" ) if $args{readtype}->{paired};
	my $READS_PIPE = ps_open( "+>$reads_file" );

	# parent process streams out sequences to fifos
	# child processes run the search on incoming sequences
	debug "Demuxing sequences\n";
	demux_sequences(
					 bowtie2_pipe1 => $BOWTIE2_R1_PIPE,
					 bowtie2_pipe2 => $BOWTIE2_R2_PIPE,
					 lastal_pipes   => \@LAST_PIPE_ARRAY,
					 dna           => $self->{"dna"},
					 file1         => $self->{"readsFile"},
					 file2         => $self->{"readsFile_2"},
					 reads_pipe    => $READS_PIPE,
	);

	# join with children when the searches are done
	foreach (@children) {
		my $tmp = waitpid( $_, 0 );
	}

	# clean up
	`rm -f $bowtie2_r1_pipe $reads_file`;
	`rm -f $bowtie2_r2_pipe` if defined($bowtie2_r2_pipe);
	foreach my $last_pipe(@last_pipe_array){
		`rm -f $last_pipe`;
	}
}

=head2 demux_sequences

reads a sequence file and streams it out to named pipes

=cut

sub demux_sequences {
	my %args           = @_;
	my $BOWTIE2_PIPE1  = $args{bowtie2_pipe1} || miss("bowtie2_pipe1");
	my $BOWTIE2_PIPE2  = $args{bowtie2_pipe2};
	my $READS_PIPE     = $args{reads_pipe} || miss("reads_pipe");
	my $last_array_reference     = $args{lastal_pipes} || miss("lastal_pipes");
	my @LAST_PIPE_ARRAY = @{$last_array_reference};
	my $F1IN           = Phylosift::Utilities::open_sequence_file( file => $args{file1} );
	my $F2IN;
	$F2IN = Phylosift::Utilities::open_sequence_file( file => $args{file2} ) if length( $args{file2} ) > 0;
	my @lines1;
	my @lines2;
	$lines1[0] = <$F1IN>;
	$lines2[0] = <$F2IN> if defined($F2IN);
	my $lastal_index=0;
	my $lastal_threads = scalar(@LAST_PIPE_ARRAY);
	while ( defined( $lines1[0] ) ) {
		my $last_pipe = $LAST_PIPE_ARRAY[$lastal_index];
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
			print $last_pipe $lines1[0] . $lines1[1] ;
			print $last_pipe $lines2[0] . $lines2[1] if defined($F2IN);
			
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

			if ( length( $lines1[1] ) > 2000 || ( defined($F2IN) && length( $lines2[1] ) > 2000 ) ) {
				# send the reads to bowtie
				print $BOWTIE2_PIPE1 @lines1;
				print $BOWTIE2_PIPE2 @lines2 if defined($F2IN);
			}

			# send the reads to last
			print $last_pipe $lines1[0] . $lines1[1] ;
			print $last_pipe $lines2[0] . $lines2[1] if defined($F2IN);

			#
			# send the reads to the reads file to write candidates later
			print $READS_PIPE $lines1[0] . $lines1[1];
			print $READS_PIPE $lines2[0] . $lines2[1] if defined($F2IN);
			@lines1 = ( $newline1, "" );
			@lines2 = ( $newline2, "" );
		}
		$lastal_index++;
		$lastal_index = $lastal_index % $lastal_threads;
	}
	
	foreach my $LAST_PIPE (@LAST_PIPE_ARRAY){
		close($LAST_PIPE);
	}
	close($BOWTIE2_PIPE1);
	close($BOWTIE2_PIPE2) if defined($F2IN);
	close($READS_PIPE);
}

sub read_marker_lengths {
	my %args = @_;
	my $self = $args{self};
	foreach my $marker (@markers) {
		$markerLength{$marker} = Phylosift::Utilities::get_marker_length( self => $self, marker => $marker );
	}
}

=head2 lastal_table

runs lastal with out-of-frame (OOF) detection on a query file (or named pipe)
returns a stream file handle

=cut

sub lastal_table {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $query_file = $args{query_file} || miss("query_file");
	my $lastal_cmd = "$Phylosift::Utilities::lastal -F15 -e75 -f0 $Phylosift::Utilities::marker_dir/replast $query_file |";
	debug "Running $lastal_cmd";
	my $HISTREAM = ps_open( $lastal_cmd );
	return $HISTREAM;
}

=head2 lastal_table_rna

runs lastal rna search on a query file (or named pipe)
returns a stream file handle

=cut

sub lastal_table_rna {
	my $self       = shift || miss("self");
	my $query_file = shift || miss("query_file");
	my $lastal_cmd = "$Phylosift::Utilities::lastal -e300 -f0 $Phylosift::Utilities::marker_dir/rnadb $query_file |";
	debug "Running $lastal_cmd";
	my $HISTREAM = ps_open( $lastal_cmd );
	return $HISTREAM;
}

=head2 bowtie2

runs bowtie2 on a single or pair of query files (or named pipes)
returns a stream file handle

=cut

sub bowtie2 {
	my %args     = @_;
	my $self     = $args{self} || miss("self");
	my $readtype = $args{readtype} || miss("readtype");
	debug "INSIDE bowtie2\n";
	my $bowtie2_cmd =
	    "$Phylosift::Utilities::bowtie2align -x "
	  . Phylosift::Utilities::get_bowtie2_db( self => $self )
	  . " --quiet --sam-nohead --sam-nosq --maxins 1000 --local ";
	$bowtie2_cmd .= " -f " if $args{readtype}->{format} eq "fasta";
	if ( $args{readtype}->{paired} ) {
		$bowtie2_cmd .= " -1 $args{reads1} -2 $args{reads2} ";
	} else {
		$bowtie2_cmd .= " -U $args{reads1} ";
	}
	$bowtie2_cmd .= " |";
	debug "Running $bowtie2_cmd";
	my $HISTREAM = ps_open( $bowtie2_cmd );
	return $HISTREAM;
}

=head2 translate_frame

=cut

sub translate_frame {
	my %args              = @_;
	my $id                = $args{id} || miss("id");
	my $seq               = $args{seq} || miss("seq");
	my $frame             = $args{frame} || miss("frame");
	my $marker            = $args{marker} || miss("marker");
	my $reverse_translate = $args{reverse_translate} || miss("reverse_translate");
	my $return_seq        = "";
	my $new_seq           = Bio::LocatableSeq->new( -seq => $seq, -id => 'temp', -verbose => 0 );
	$new_seq = $new_seq->revcom() if ( $frame < 0 );

	if ($reverse_translate) {
		$id = Phylosift::Summarize::tree_name( name => $id );
		if ( exists $markerNuc{$marker} ) {
			$markerNuc{$marker} .= ">" . $id . "\n" . $new_seq->seq . "\n";
		} else {
			$markerNuc{$marker} = ">" . $id . "\n" . $new_seq->seq . "\n";
		}
	}
	$return_seq = $new_seq->translate();
	return $return_seq->seq();
}

=head2 get_hits_contigs

parse the blast file

=cut

sub get_hits_contigs {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $HITSTREAM  = $args{HITSTREAM} || miss("HITSTREAM");
	my $searchtype = $args{searchtype} || miss("searchtype");    # can be blastx or lastal

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
		my ( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore, $frameshift );
		if ( $searchtype eq "blastx" ) {
			( $query, $subject, $two, $three, $four, $five, $query_start, $query_end, $eight, $nine, $ten, $bitScore ) = split( /\t/, $_ );
		} else {

			# read table in lastal format
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
			$frameshift = $dat[11];
		}

		# get the marker name
		my @marker = split( /\_\_/, $subject );    # this is soooo ugly
		my $markerName = $marker[0];

		# running on long reads or an assembly
		# allow each region of a sequence to have a top hit
		# do not allow overlap
		if ( defined( $contig_top_bitscore{$query} ) ) {
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

					# debug "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end, $frameshift ] if ( $bitScore > $prevhit[1] );
					last;
				}

				# now check the same for reverse-strand hits
				if (    $prevhit[2] > $prevhit[3]
					 && $query_start > $query_end
					 && $prevhit[3] < $query_start - $max_hit_overlap
					 && $query_end + $max_hit_overlap < $prevhit[2] )
				{

					# debug "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end, $frameshift ] if ( $bitScore > $prevhit[1] );
					last;
				}
			}
			if ( $i == @{ $contig_hits{$query} } ) {

				# no overlap was found, include this hit
				my @hitdata = [ $markerName, $bitScore, $query_start, $query_end, $frameshift ];
				push( @{ $contig_hits{$query} }, @hitdata );
			}
		} elsif ( !defined( $contig_top_bitscore{$query} ) ) {
			my @hitdata = [ $markerName, $bitScore, $query_start, $query_end, $frameshift ];
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
	my $self       = $args{self} || miss("self");
	my $HITSTREAM  = $args{HITSTREAM} || miss("HISTREAM");
	my $searchtype = $args{searchtype} || miss("searchtype");
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
		if ( $self->{"isolate"} ) {

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
	my $self      = $args{self} || miss("self");
	my $HITSTREAM = $args{HITSTREAM} || miss("HISTREAM");
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
	my $subject     = $args{subject} || miss("subject");
	my $search_type = $args{search_type} || miss("search_type");
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

=head2 write_candidates

write out results

=cut

sub write_candidates {
	my %args          = @_;
	my $self          = $args{self} || miss("self");
	my $contigHitsRef = $args{hitsref} || miss("hitsref");
	my $type       = $args{searchtype} || ""; # search type -- candidate filenames will have this name embedded, enables parallel output from different programs
	my $reads_file = $args{reads} || miss("reads");
	my $process_id = $args{process_id} || miss("process id");
	my %contig_hits = %$contigHitsRef;
	my %markerHits;
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

			# ensure flanking region is a multiple of 3 to avoid breaking frame in DNA
			my $seqLength = length( $seq->seq );
			my $newSeq;
			$newSeq = substr( $seq->seq, $start, $end - $start );

			#if we're working from DNA then need to translate to protein
			if ( $self->{"dna"} && $type !~ /\.rna/ ) {
				# check for frameshifts
				# frameshift from lastal is a string of <block>,<gaps>,<block>
				# where blocks are aligned chunks and gaps are of the form X:Y
				# listing number of gaps in ref and query seq, respectively
				my $frameshift = $cur_hit[4];
#				debug $seq->id ." Frameshift $frameshift, start $start, end $end, seqlen ".$seq->length."\n";
				my @frames = split(/,/, $frameshift);
				
				# compute alignment length across all frames
				my $aln_len = 0;
				foreach my $frame(@frames){
					next if $frame =~ /:/;	# ignore gap regions
					$aln_len += $frame;
				}
				# check length again in AA units
				$min_len = $markerLength{$markerHit} < $seq->length / 3 ? $markerLength{$markerHit} : $seq->length / 3;
				next unless ( $aln_len / $min_len >= $align_fraction );
				
				# construct the DNA sequence to translate
				my $dna_seq = "";
				my $cs = $start;
				my $fstart = $cur_hit[2] > $cur_hit[3] ? scalar(@frames)-1 : 0;
				my $fdiff = $cur_hit[2] > $cur_hit[3] ? -1 : 1;
				for(my $fI=$fstart; $fI >= 0 && $fI < @frames; $fI += $fdiff){
					my $frame = $frames[$fI];
					if($frame =~ /(\d+):-*(\d+)/){
						my $fs = $2;
						$fs *= -1 if $frame =~ /:-/;
						$cs += $fs;
#						debug "Found frameshift, $fs bases\n";
					}else{
						$dna_seq .= $seq->subseq($cs, $cs+($frame*3)-1);
						$cs += $frame*3;
					}
				}

				my $strand = 1;	# forward or reverse strand?
				$strand *= -1 if ( $cur_hit[2] > $cur_hit[3] );

				$newSeq = translate_frame(
										   id                => $seq->id,
										   seq               => $dna_seq,
										   frame             => $strand,
										   marker            => $markerHit,
										   reverse_translate => $self->{"dna"}
				);
				$newSeq =~ s/\*/X/g;    # bioperl uses * for stop codons but we want to give X to hmmer later
			}
			$markerHits{$markerHit} = "" unless defined( $markerHits{$markerHit} );
			$markerHits{$markerHit} .= ">" . $seq->id . "\n" . $newSeq . "\n";
		}
	}

	#write the read+ref_seqs for each markers in the list
	foreach my $marker ( keys %markerHits ) {

		#writing the hits to the candidate file
		my $candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type );
		my $FILE_OUT = ps_open( ">$candidate_file".".".$process_id );
		print $FILE_OUT $markerHits{$marker};
		close($FILE_OUT);
		if ( $self->{"dna"} && $type !~ /\.rna/ ) {
			$candidate_file = Phylosift::Utilities::get_candidate_file( self => $self, marker => $marker, type => $type, dna => 1 );
			my $FILE_OUT = ps_open(">$candidate_file".".".$process_id );
			print $FILE_OUT $markerNuc{$marker} if defined( $markerNuc{$marker} );
			close($FILE_OUT);
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
	my $self = $args{self} || miss("self");
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
