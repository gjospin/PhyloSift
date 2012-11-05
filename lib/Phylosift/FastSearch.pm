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
use Phylosift::MarkerAlign;
use Phylosift::Settings;
use File::Basename;
use POSIX qw(ceil floor);

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

set_default_values();
my $best_hits_bit_score_range = $Phylosift::Settings::best_hits_bit_score_range;

# all hits with a bit score within this amount of the best will be used
my $align_fraction = $Phylosift::Settings::align_fraction;

# at least this amount of min[length(query),length(marker)] must align to be considered a hit
my $align_fraction_isolate = $Phylosift::Settings::align_fraction_isolate;

# use this align_fraction when in isolate mode on long sequences
my $quality_threshold = $Phylosift::Settings::quality_threshold;
my @lookup_array      = ();
my %markers;
my %markerNuc = ();
my %markerLength;

sub run_search {
	my %args       = @_;
	my $self       = $args{self} || miss("PS object");
	my $markersRef = $args{marker_reference} || miss("marker_reference");
	%markers = map { $_ => 1 } @{$markersRef};
	debug "USING ".keys(%markers)."\n";
	my $position = rindex( $self->{"readsFile"}, "/" );
	$self->{"readsFile"} =~ m/(\w+)\.?(\w*)$/;

	#setting default values for the module
	#set_default_values(self=>$self);

	# set align_fraction appropriately
	$align_fraction = $align_fraction_isolate if ($Phylosift::Settings::isolate);

	my $input_file = $self->{stdin} ? "STDIN" : $self->{"readsFile"};
	my $F1IN = Phylosift::Utilities::open_sequence_file( file => $input_file );

	# check what kind of input was provided
	my $type = Phylosift::Utilities::get_sequence_input_type($F1IN);
	$self->{readtype} = $type;
	$self->{"dna"} = $type->{seqtype} eq "protein" ? 0 : 1;    # Is the input protein sequences?
	debug "Input type is $type->{seqtype}, $type->{format}\n";

	#making sure $type->{paired} is set so we create the appropriate variables
	$type->{paired} = 1 if ( exists $self->{"readsFile_2"} && length( $self->{"readsFile_2"} ) > 0 );

	# ensure databases and sequences are prepared for search
	prep_and_clean( self => $self );

	#read_marker_lengths( self => $self );
	# search reads/contigs against marker database
	my $searchtype = "blast";
	my $contigs    = 0;
	$contigs = 1
	  if (    defined($Phylosift::Settings::coverage)
		   && $type->{seqtype} ne "protein"
		   && ( !defined $Phylosift::Settings::isolate || $Phylosift::Settings::isolate != 1 ) );

	# open the input file(s)
	my $F2IN;
	$F2IN = Phylosift::Utilities::open_sequence_file( file => $self->{"readsFile_2"} )
	  if length( $self->{"readsFile_2"} ) > 0;
	my $RUNINFO = ps_open( ">>".Phylosift::Utilities::get_run_info_file( self => $self ) );

	# launch the searches
	my $start_chunk = defined($Phylosift::Settings::start_chunk) ? $Phylosift::Settings::start_chunk : 1;
	for ( my $chunkI = 1;; $chunkI++ ) {

		#reads the marker_summary.txt from blastDir

		my $completed_chunk = has_chunk_completed( self => $self, chunk => $chunkI );
		$completed_chunk = 1 if $start_chunk > $chunkI;

		# need to run this even if chunk is done so that we advance through the input file
		# TODO: don't launch lastal on RNA unless really needed
		my $finished = launch_searches(
										self     => $self,
										readtype => $type,
										dir      => $self->{"blastDir"},
										contigs  => $contigs,
										chunk    => $chunkI,
										FILE1    => $F1IN,
										FILE2    => $F2IN
		);
		compute_hits_summary( self => $self, chunk => $chunkI );
		if ( !$completed_chunk && ( $self->{"mode"} eq "all" || $self->{"continue"} ) ) {

			# fire up the next step!
			# TODO: make this a call to a function "chunk_done" in main Phylosift module
			# that starts the next step
			Phylosift::Utilities::end_timer( name => "runBlast" );
			Phylosift::Utilities::start_timer( name => "runAlign" );
			Phylosift::MarkerAlign::MarkerAlign(
												 self             => $self,
												 marker_reference => $markersRef,
												 chunk            => $chunkI
			);
		}
		print $RUNINFO "Chunk $chunkI completed\n" unless $completed_chunk;
		debug "Debug lvl : $Phylosift::Utilities::debuglevel\n";
		clean_chunk_directory( self => $self, chunk => $chunkI ) if !$Phylosift::Settings::keep_search;
		last if $finished || ( defined($Phylosift::Settings::chunks) && ( $chunkI - $start_chunk + 1 ) >= $Phylosift::Settings::chunks );
	}
	return $self;
}

=head2 compute_hits_summary

	reads all candidate files and compiles hit numbers for each marker

=cut

sub compute_hits_summary {
	my %args                    = @_;
	my $chunk                   = $args{chunk} || miss("chunk");
	my $self                    = $args{self} || miss("PS Object");
	my $marker_hits_numbers_ref = Phylosift::Utilities::read_marker_summary( self => $self, path => $self->{"blastDir"} );
	my %marker_hits_numbers     = %{$marker_hits_numbers_ref};
	my @candidates              = glob( $self->{"blastDir"}."/*.aa.$chunk*" );
	foreach my $cand_file (@candidates) {
		my $grep_cmd = "grep '>' -c $cand_file";
		$cand_file =~ m/blastDir\/([^\.]+)/;
		my $grep_results = `$grep_cmd`;
		chomp($grep_results);
		$marker_hits_numbers{$1} = 0 unless exists $marker_hits_numbers{$1};
		$marker_hits_numbers{$1} += $grep_results;
	}
	Phylosift::Utilities::print_marker_summary( self => $self, path => $self->{"blastDir"}, summary => \%marker_hits_numbers );
}

=head2 set_default_values

	Set the necessary default values for this module if they haven't been read in by the RC file.

=cut

sub set_default_values {
	my %args               = @_;
	my $default_chunk_size = 200000;
	$default_chunk_size = 1000000 if !defined($Phylosift::Settings::extended) || !$Phylosift::Settings::extended;
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::CHUNK_MAX_SEQS,            value => $default_chunk_size );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::lastal_evalue,             value => "-e75" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::lastal_rna_evalue,         value => "-e75" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::max_hit_overlap,           value => 10 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::discard_length,            value => 30 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::best_hits_bit_score_range, value => 30 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::align_fraction,            value => 0.5 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::align_fraction_isolate,    value => 0.8 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::quality_threshold,         value => 20 );
}

=head2 clean_search_directory
	
	Deleting the files from the blastDir for a specific chunk
	Passing the chunk so this is ready in case PS moves towards running multiples chunks in parallel.

=cut

sub clean_chunk_directory {
	my %args  = @_;
	my $self  = $args{self};
	my $chunk = $args{chunk};
	return if $Phylosift::Utilities::debuglevel >= 1;    #Do not clean if debug is present
	my $remove_aa  = "rm ".$self->{"blastDir"}."/*.aa.$chunk*";
	my $remove_ffn = "rm ".$self->{"blastDir"}."/*.ffn.$chunk*";
	debug("Cleaning the Search directory for chunk $chunk\n");
	my @array_to_delete = glob( $self->{"blastDir"}."/*.ffn.$chunk*" );
	`$remove_aa`;                                        #deletes the aa files
	`$remove_ffn` if @array_to_delete;                   #added check to prevent an error when using AA sequences;
}

=head2 has_chunk_completed

Checks to see if a chunk has completed
returns 1 if it has and 0 if it hasn't
=cut

sub has_chunk_completed {
	my %args     = @_;
	my $self     = $args{self} || miss("PS object");
	my $chunk    = $args{chunk};
	my $run_file = Phylosift::Utilities::get_run_info_file( self => $self );
	my $grep     = `grep "Chunk $chunk completed" $run_file`;
	return 1 if defined $grep && length($grep) > 0 && !$Phylosift::Settings::force;
	return 0;
}

=head2 launch_searches

creates named pipes to stream input to search programs and launches them
after creating pipes, the program forks into separate processes.
the parent process writes sequence data to files, child processes launch a similarity search on that data

=cut

sub launch_searches {
	my %args          = @_;
	my $self          = $args{self};
	my $dir           = $args{dir} || miss("dir");
	my $readtype      = $args{readtype} || miss("readtype");
	my $chunk         = $args{chunk} || miss("chunk");
	my $FILE1         = $args{FILE1};
	my $FILE2         = $args{FILE2};
	my $chunky        = defined($chunk) ? ".$chunk" : "";
	my $reads_file    = $dir."/reads.fasta$chunky";
	my $last_rna_pipe = $dir."/last_rna.pipe";
	debug "Making fifos\n";
	my @last_pipe_array = ();

	for ( my $i = 0; $i <= $Phylosift::Settings::threads - 1; $i++ ) {
		push( @last_pipe_array, $dir."/last_$i.pipe" );
		`mkfifo "$dir/last_$i.pipe"`;
	}
	my $rna_procs = 0;
	if ( $readtype->{seqtype} eq "dna"
		 && -e Phylosift::Utilities::get_rna_db( self => $self ).".prj" )
	{
		`mkfifo "$last_rna_pipe"`;
		$rna_procs = 1;
	} else {

		# if we're searching protein, send these to the bitbucket
		$last_rna_pipe = "/dev/null";
	}
	my @children;
	for ( my $count = 1; $count <= $Phylosift::Settings::threads + $rna_procs; $count++ ) {
		my $pid = fork();
		if ($pid) {

			# parent process will write sequences below
			push( @children, $pid );
		} elsif ( $pid == 0 ) {
			debug "Launching search process $count\n";

			# child processes will search sequences
			my $hitstream;
			my $candidate_type = ".$count";
			if ( $count <= $Phylosift::Settings::threads ) {
				$hitstream = lastal_table(
										   self       => $self,
										   query_file => $last_pipe_array[ $count - 1 ],
										   readtype   => $readtype
				);
				$candidate_type = ".lastal";
			} elsif ( $count == $Phylosift::Settings::threads + 1 ) {

				if ( !-e Phylosift::Utilities::get_rna_db( self => $self ).".prj" ) {
					debug "Exiting process $count because rna db not found\n";
					my $lrp1 = ps_open($last_rna_pipe);
					`rm -f $last_rna_pipe`;
					exit 0;
				} else {
					$hitstream = lastal_table_rna( self       => $self,
												   query_file => $last_rna_pipe );
				}
				$candidate_type = ".lastal.rna";
			}
			my $hitsref = get_hits_contigs(
											self       => $self,
											HITSTREAM  => $hitstream,
											searchtype => "lastal",
											pid        => $count
			);

			# write out sequence regions hitting marker genes to candidate files
			debug "Writing candidates from process $count\n";
			write_candidates(
							  self       => $self,
							  hitsref    => $hitsref,
							  searchtype => "$candidate_type",
							  reads      => $reads_file,
							  process_id => $count,
							  chunk      => $chunk
			);
			exit 0;
		} else {
			croak "couldn't fork: $!\n";
		}
	}
	my @LAST_PIPE_ARRAY = ();
	foreach my $last_pipe (@last_pipe_array) {
		push( @LAST_PIPE_ARRAY, ps_open(">$last_pipe") );
	}

	debug "Opening $reads_file\n";
	my $READS_PIPE    = ps_open("+>$reads_file");
	my $LAST_RNA_PIPE = ps_open(">$last_rna_pipe");

	# parent process streams out sequences to fifos
	# child processes run the search on incoming sequences
	debug "Octopus is handing out sequences\n";
	my $finished = demux_sequences(
									lastal_pipes  => \@LAST_PIPE_ARRAY,
									LAST_RNA_PIPE => $LAST_RNA_PIPE,
									dna           => $self->{"dna"},
									FILE1         => $FILE1,
									FILE2         => $FILE2,
									reads_pipe    => $READS_PIPE,
									readtype      => $readtype,
									chunk         => $chunk,
									self          => $self,
	);

	# join with children when the searches are done
	foreach (@children) {
		my $tmp = waitpid( $_, 0 );
	}

	# clean up
	`rm -f "$reads_file"`;
	`rm -f "$last_rna_pipe"`;
	foreach my $last_pipe (@last_pipe_array) {
		`rm -f "$last_pipe"`;
	}
	return $finished;
}

=head2 get_next_line

Read a line from an array-buffered input stream. 
If there a line is available in the buffer then it is returned and the buffer index is incremented.
Otherwise the next line from the input stream is returned. 

=cut

sub get_next_line {
	my %args         = @_;
	my $buffer       = $args{buffer} || miss("buffer");
	my $buffer_index = $args{buffer_index};
	my $FILE         = $args{FILE} || miss("FILE");
	my $line;
	if ( defined( @$buffer[$$buffer_index] ) ) {
		$line = @$buffer[$$buffer_index];
		@$buffer[ $$buffer_index++ ] = undef;
	} else {
		$line = <$FILE>;
	}
	return $line;
}

=head2 demux_sequences

reads a sequence file and streams it out to named pipes

=cut

sub demux_sequences {
	my %args                 = @_;
	my $READS_PIPE           = $args{reads_pipe} || miss("reads_pipe");
	my $last_array_reference = $args{lastal_pipes} || miss("lastal_pipes");
	my $self                 = $args{self} || miss("PS object");
	my $chunk                = $args{chunk} || miss("Chunk number");
	my $LAST_RNA_PIPE        = $args{LAST_RNA_PIPE};
	my @LAST_PIPE_ARRAY      = @{$last_array_reference};
	my $F1IN                 = $args{FILE1};
	my $F2IN                 = $args{FILE2};
	my $paired               = $Phylosift::Settings::paired;
	my $readtype             = $args{readtype};
	my @lines1;
	my @lines2;
	my $lookup_id_filename = $self->{"blastDir"}."/lookup_ID.$chunk.tbl";
	my $IDFILE             = ps_open(">$lookup_id_filename");
	my $lastal_index       = 0;
	my $lastal_threads     = scalar(@LAST_PIPE_ARRAY);
	my $completed_chunk    = has_chunk_completed( self => $self, chunk => $chunk );

	# the following two variables track how far we are through a chunk
	my $seq_count    = 0;
	my $seq_size     = 0;
	my $buffer_index = 0;
	if ( defined( $self->{"stashed_lines"} ) ) {
		$lines1[0] = $self->{"stashed_lines"}[0];
		$lines2[0] = $self->{"stashed_lines"}[1] if defined( $self->{"stashed_lines"}[1] );
	} else {
		$lines1[0] = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN );
		$lines2[0] = <$F2IN> if defined($F2IN);
	}

	while ( defined( $lines1[0] ) ) {
		$seq_count++;
		my $last_pipe = $LAST_PIPE_ARRAY[$lastal_index];
		if ( $lines1[0] =~ /^@/ ) {

			#
			# FASTQ format
			for ( my $i = 1; $i < 4; $i++ ) {
				$lines1[$i] = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN );
				$lines2[$i] = <$F2IN> if defined($F2IN);
			}
			if ( $paired && !defined($F2IN) ) {
				for ( my $i = 0; $i < 4; $i++ ) {
					$lines2[$i] = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN );
				}
			}

			# quality trim the read(s)
			qtrim_read( read => \@lines1, quality => $quality_threshold, readtype => $readtype );
			qtrim_read( read => \@lines2, quality => $quality_threshold, readtype => $readtype ) if defined( $lines2[0] );

			#add the reads to file lookup

			#			if ( $lines1[0] =~ m/^@(\S+)(\/\d)/ && $paired ) {
			chomp( $lines1[0] );
			print $IDFILE "$lines1[0]\t$seq_count/1\n";

			#			} elsif ( $lines1[0] =~ m/^@(.+)/ ) {
			if ( defined( $lines2[0] ) ) {
				chomp( $lines2[0] );
				print $IDFILE "$lines2[0]\t$seq_count/2\n";
			}

			$lines1[0] = "\@$seq_count/1\n";
			if ($paired) {
				$lines2[0] = "\@$seq_count";
				$lines2[0] .= "/2" if $paired;
				$lines2[0] .= "\n";
			}

			# some FastQ use . instead of N (wtf?)
			$lines1[1] =~ s/\./N/g;
			$lines2[1] =~ s/\./N/g if @lines2;

			print $LAST_RNA_PIPE @lines1 unless $completed_chunk;
			print $LAST_RNA_PIPE @lines2 if @lines2 && !$completed_chunk;

			#
			# send the reads to lastal (convert to fasta)
			$lines1[0] =~ s/^@/>/g;
			$lines2[0] =~ s/^@/>/g if @lines2;

			print $last_pipe $lines1[0].$lines1[1] unless $completed_chunk;
			print $last_pipe $lines2[0].$lines2[1] if @lines2 && !$completed_chunk;

			#
			# send the reads to the reads file in fasta format to write candidates later
			print $READS_PIPE $lines1[0].$lines1[1] unless $completed_chunk;
			print $READS_PIPE $lines2[0].$lines2[1] if @lines2 && !$completed_chunk;

			#
			# prepare for next loop iter
			$lines1[0] = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN );
			$lines2[0] = <$F2IN> if defined($F2IN);
		} elsif ( $lines1[0] =~ /^>/ ) {

			#
			# FASTA format
			my $newline1;
			while ( $newline1 = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN ) ) {
				last if $newline1 =~ /^>/;
				chomp( $lines1[1] ) if defined( $lines1[1] );
				$lines1[1] .= $newline1;
			}
			my $newline2;
			if ( defined($F2IN) ) {
				while ( $newline2 = <$F2IN> ) {
					last if $newline2 =~ /^>/;
					chomp( $lines2[1] ) if defined( $lines2[1] );
					$lines2[1] .= $newline2;
				}
			}
			if ( $paired && !defined($F2IN) ) {
				@lines2 = ( $newline1, "" );
				while ( $newline1 = get_next_line( buffer => $readtype->{buffer}, buffer_index => \$buffer_index, FILE => $F1IN ) ) {
					last if $newline1 =~ /^>/;
					chomp( $lines1[1] ) if defined( $lines1[1] );
					$lines2[1] .= $newline1;
				}
			}

			#removing possible \n in the middle of the sequence
			$lines1[1] =~ s/\n//g;
			$lines2[1] =~ s/\n//g if defined( $lines2[1] );

			#add a \n at the end of the sequence to comply with fasta format
			$lines1[1] .= "\n";
			$lines2[1] .= "\n" if defined( $lines2[1] );

			#add the reads to file lookup
			#			if ( $lines1[0] =~ m/^>(\S+)(\/\d)/ && $paired ) {
			$lines1[0] =~ s/^>//;
			chomp( $lines1[0] );
			print $IDFILE "$lines1[0]\t$seq_count/1\n";

			#			} elsif ( $lines1[0] =~ m/^>(.+)/ ) {
			if ( defined $lines2[0] ) {
				$lines2[0] =~ s/^>//;
				chomp( $lines2[0] );
				print $IDFILE "$lines2[0]\t$seq_count/2\n" if defined( $lines2[0] );
			}

			$lines1[0] = ">$seq_count/1\n";
			if ($paired) {
				$lines2[0] = ">$seq_count";
				$lines2[0] .= "/2" if $paired;
				$lines2[0] .= "\n";
			}

			# if reads are long, do RNA search with lastal
			print $LAST_RNA_PIPE $lines1[0].$lines1[1] unless $completed_chunk;
			print $LAST_RNA_PIPE $lines2[0].$lines2[1] if @lines2 && !$completed_chunk;

			# send the reads to last
			print $last_pipe $lines1[0].$lines1[1] unless $completed_chunk;
			print $last_pipe $lines2[0].$lines2[1] if @lines2 && !$completed_chunk;

			#
			# send the reads to the reads file to write candidates later
			print $READS_PIPE $lines1[0].$lines1[1] unless $completed_chunk;
			print $READS_PIPE $lines2[0].$lines2[1] if @lines2 && !$completed_chunk;

			@lines1 = ( $newline1, "" );
			@lines2 = ( $newline2, "" ) if @lines2 && !$completed_chunk;
		}
		$lastal_index++;
		$lastal_index = $lastal_index % $lastal_threads;
		last if ( $seq_count >= $Phylosift::Settings::CHUNK_MAX_SEQS );
	}
	foreach my $LAST_PIPE (@LAST_PIPE_ARRAY) {
		close($LAST_PIPE);
	}
	close($IDFILE);

	#	close($BOWTIE2_PIPE2) if defined($F2IN);
	close($LAST_RNA_PIPE);
	close($READS_PIPE);
	debug "Octopus handed out $seq_count sequences\n";

	# if there is more sequence, save the header we currently have for later
	if ( defined( $lines1[0] ) ) {
		$self->{"stashed_lines"}[0] = $lines1[0];
		$self->{"stashed_lines"}[1] = $lines2[0] if defined( $lines2[0] );
	}

	# return 1 if we've processed all the data, zero otherwise
	return 0 if defined( $lines1[0] );
	return 1;
}

=head2 qtrim_read

trims a fastq read to a particular quality score using Heng Li's algorithm from bwa.
code based on SGA's implementation.

=cut

sub qtrim_read {
	my %args     = @_;
	my $read     = $args{read};
	my $q        = $args{quality};
	my $readtype = $args{readtype};

	$q += 33 if $readtype->{qtype} eq "phred33";
	$q += 64 if $readtype->{qtype} eq "phred64";

	# Perform a soft-clipping of the sequence by removing low quality bases from the
	# 3' end using Heng Li's algorithm from bwa

	my $seq = @$read[1];
	chomp $seq;
	my $qq = @$read[3];
	chomp $qq;
	my @qual          = split( //, $qq );
	my $endpoint      = 0;                  # not inclusive
	my $max           = 0;
	my $i             = length($seq) - 1;
	my $terminalScore = ord( $qual[$i] );

	# Only perform soft-clipping if the last base has qual less than $q
	return if ( $terminalScore >= $q );

	my $subSum = 0;
	while ( $i >= 0 ) {
		my $ps    = ord( $qual[$i] );
		my $score = $q - $ps;
		$subSum += $score;
		if ( $subSum > $max ) {
			$max      = $subSum;
			$endpoint = $i;
		}
		$i--;
	}

	# Clip the read
	@$read[1] = substr( $seq,      0, $endpoint )."\n";
	@$read[3] = substr( @$read[3], 0, $endpoint )."\n";
}

sub read_marker_lengths {
	my %args = @_;
	my $self = $args{self};
	foreach my $marker ( keys(%markers) ) {
		$markerLength{$marker} = Phylosift::Utilities::get_marker_length( self   => $self,
																		  marker => $marker );
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
	my $readtype   = $args{readtype} || miss("readtype");
	my $last_opts  = "";
	$last_opts = "-F15" if $readtype->{seqtype} eq "dna";
	my $db = Phylosift::Utilities::get_lastal_db( self => $self );

	#my $lastal_cmd =
	#  "$Phylosift::Settings::lastal $last_opts $Phylosift::Settings::lastal_evalue -f0 $db \"$query_file\" |";

	my $lastal_cmd = "$Phylosift::Settings::lastal $last_opts $Phylosift::Settings::lastal_evalue ";

	#$lastal_cmd .= "-Q2 " if $self->{readtype}{format} eq "fastq";
	$lastal_cmd .= "-f0 $db \"$query_file\" |";
	debug "Running $lastal_cmd";
	my $HISTREAM = ps_open($lastal_cmd);
	return $HISTREAM;
}

=head2 lastal_table_rna

runs lastal rna search on a query file (or named pipe)
returns a stream file handle

=cut

sub lastal_table_rna {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $query_file = $args{query_file} || miss("query_file");
	my $db         = Phylosift::Utilities::get_rna_db( self => $self );

	#my $lastal_cmd =
	#"$Phylosift::Settings::lastal $Phylosift::Settings::lastal_rna_evalue -f0 $db \"$query_file\" |";
	my $lastal_cmd = "$Phylosift::Settings::lastal $Phylosift::Settings::lastal_rna_evalue ";
	$lastal_cmd .= "-Q2 " if $self->{readtype}{format} eq "fastq";
	$lastal_cmd .= "-f0 $db \"$query_file\" |";
	debug "Running $lastal_cmd";
	my $HISTREAM = ps_open($lastal_cmd);
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
	my $reverse_translate = $args{reverse_translate}
	  || miss("reverse_translate");
	my $return_seq = "";
	my $new_seq = Bio::LocatableSeq->new( -seq => $seq, -id => 'temp', -verbose => 0 );
	$new_seq = $new_seq->revcom() if ( $frame < 0 );

	if ($reverse_translate) {

		#$id = Phylosift::Summarize::tree_name( name => $id );
		if ( exists $markerNuc{$marker} ) {
			$markerNuc{$marker} .= ">".$id."\n".$new_seq->seq."\n";
		} else {
			$markerNuc{$marker} = ">".$id."\n".$new_seq->seq."\n";
		}
	}
	$return_seq = $new_seq->translate();
	return $return_seq->seq();
}

sub is_overlapping {
	my %args            = @_;
	my $prevhit         = $args{prevhit} || miss("prevhit");
	my $start           = $args{start} || miss("start");
	my $end             = $args{end} || miss("end");
	my $max_hit_overlap = $Phylosift::Settings::max_hit_overlap;

	my $s   = $start < $end                  ? $start        : $end;
	my $e   = $start >= $end                 ? $start        : $end;
	my $ph2 = $prevhit->[2] < $prevhit->[3]  ? $prevhit->[2] : $prevhit->[3];
	my $ph3 = $prevhit->[2] >= $prevhit->[3] ? $prevhit->[2] : $prevhit->[3];

	return 1 if ( $ph2 < $e - $max_hit_overlap && $s + $max_hit_overlap < $ph3 );
	return 0;
}

=head2 get_hits_contigs

parse the LASTAL results

=cut

sub get_hits_contigs {
	my %args             = @_;
	my $self             = $args{self} || miss("self");
	my $HITSTREAM        = $args{HITSTREAM} || miss("HITSTREAM");
	my $searchtype       = $args{searchtype} || miss("searchtype");    # can be blastx or lastal
	my $pid              = $args{pid} || miss("pid");
	my $prev_hit         = "";
	my %hit_counts       = ();
	my $hit_counts_total = 0;

	# key is a contig name
	# value is an array of arrays, each one has [marker,bit_score,left-end,right-end,suffix]
	my %contig_hits;
	my %contig_top_bitscore;

	# return empty if there is no data
	return \%contig_hits unless defined( fileno $HITSTREAM );
	while (<$HITSTREAM>) {

		# read a lastal line
		next if ( $_ =~ /^#/ );
		chomp($_);
		my ( $query, $subject, $query_start, $query_end, $bitScore, $frameshift );

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

		# get the marker name
		my @marker = split( /\_\_/, $subject );    # this is soooo ugly
		my $markerName = $marker[0];
		next unless ( exists $markers{$markerName} );

		# running on long reads or an assembly
		# allow each region of a sequence to have a top hit
		# do not allow overlap

		my $i   = 0;
		my $add = 0;
		if ( defined( $contig_top_bitscore{$query} ) ) {
			my $cur_hitcount = scalar( @{ $contig_hits{$query} } );
			for ( ; $i < $cur_hitcount; $i++ ) {
				my $prevhitref = $contig_hits{$query}->[$i];

				# is there enough overlap to consider these the same?
				# if so, take the new one if it has higher bitscore
				if ( is_overlapping( prevhit => $prevhitref, start => $query_start, end => $query_end ) ) {
					$add = 1 if $bitScore > $prevhitref->[1];
					last;
				}
			}

			# add if we didn't find any overlap
			$add = 1 if ( $i == $cur_hitcount );
		} elsif ( !defined( $contig_top_bitscore{$query} ) ) {
			$add = 1;
			$contig_top_bitscore{$query}{$markerName} = $bitScore;
		}
		if ($add) {

			# include this hit
			$contig_hits{$query}->[$i] = [ $markerName, $bitScore, $query_start, $query_end, $frameshift, $pid ];
			$hit_counts{$markerName} = 0 unless exists $hit_counts{$markerName};
			$hit_counts{$markerName}++;
			$hit_counts_total++;
		}
	}
	foreach my $marker ( keys %hit_counts ) {
		if ( $hit_counts_total > $Phylosift::Settings::CHUNK_MAX_SEQS / 10 && $hit_counts{$marker} / $hit_counts_total > 0.5 ) {
			if ( $self->{"custom_chunk_size"} == 0 ) {
				my $warning = "\n\n\nDetected large number of hits to the marker $marker with Default chunk size used.\n";
				$warning .= "This will require a large amount of memory. \n";
				$warning .= "Manually tune the number of sequences to be processed per chunks in phylosiftrc\n\n\n";
				Phylosift::Utilities::print_bat_signal();
				croak $warning;
			}
		}
	}
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

	# search type -- candidate filenames will have this name embedded, enables parallel output from different programs
	my $type       = $args{searchtype} || "";
	my $reads_file = $args{reads}      || miss("reads");
	my $process_id = $args{process_id} || miss("process id");
	my $chunk      = $args{chunk}      || miss("chunk");
	my %contig_hits = %$contigHitsRef;
	my %markerHits;
	my $align_fraction = $Phylosift::Settings::align_fraction;
	debug "ReadsFile:  $self->{\"readsFile\"}"."\n";
	my $SEQIN = ps_open($reads_file);

	# input file is FastA with one sequence per line
	# parsing is very simple, just alternate lines.
	while ( my $idline = <$SEQIN> ) {
		chomp $idline;
		my $id = $1 if $idline =~ /^>(.+)/;
		my $seq = <$SEQIN>;
		chomp $seq;

		# skip this one if there are no hits
		next unless ( exists $contig_hits{$id} );
		my $cur_hit_count = scalar( @{ $contig_hits{$id} } );
		for ( my $i = 0; $i < $cur_hit_count; $i++ ) {
			my $cur_hit   = $contig_hits{$id}->[$i];
			my $markerHit = $cur_hit->[0];
			my $start     = $cur_hit->[2];
			my $end       = $cur_hit->[3];
			my $suff      = $cur_hit->[5];
			( $start, $end ) = ( $end, $start ) if ( $start > $end );    # swap if start bigger than end

			# check to ensure hit covers enough of the marker
			# TODO: make this smarter about boundaries, e.g. allow a smaller fraction to hit
			# if it looks like the query seq goes off the marker boundary
			$markerLength{$markerHit} = Phylosift::Utilities::get_marker_length( self => $self, marker => $markerHit ) unless exists $markerLength{$markerHit};
			if ( !defined($markerHit) || !defined( $markerLength{$markerHit} ) ) {
				debug "no alignment length for marker $markerHit\n";
				next;
			}
			my $min_len = $markerLength{$markerHit} < length($seq) ? $markerLength{$markerHit} : length($seq);
			next unless ( ( $end - $start ) / $min_len >= $align_fraction );

			# ensure flanking region is a multiple of 3 to avoid breaking frame in DNA
			my $seqLength = length($seq);
			my $new_seq;
			$new_seq = substr( $seq, $start, $end - $start );
			my $new_id = $id;
			my $coord  = ".$start.$end";

			#			if($Phylosift::Settings::paired){
			#				#when working with paired data, move the mateID to the end of the string
			#				$new_id =~ s/(\d+)(\/\d)/$1$coord$2/;
			#			} else {
			$new_id =~ s/^(\d+)\/(\d)/$1.$2$coord/;

			#			}

			#if we're working from DNA then need to translate to protein
			if ( $self->{"dna"} && $type !~ /\.rna/ ) {

				# check for frameshifts
				# frameshift from lastal is a string of <block>,<gaps>,<block>
				# where blocks are aligned chunks and gaps are of the form X:Y
				# listing number of gaps in ref and query seq, respectively
				my $frameshift = $cur_hit->[4];

				# debug $seq->id ." Frameshift $frameshift, start $start, end $end, seqlen ".$seq->length."\n";
				my @frames = split( /,/, $frameshift );

				# compute alignment length across all frames
				my $aln_len = 0;
				foreach my $frame (@frames) {
					next if $frame =~ /:/;    # ignore gap regions
					$aln_len += $frame;
				}

				# check length again in AA units
				$min_len =
				    $markerLength{$markerHit} < length($seq) / 3
				  ? $markerLength{$markerHit}
				  : length($seq) / 3;
				next unless ( $aln_len / $min_len >= $align_fraction );

				# construct the DNA sequence to translate
				my $dna_seq = "";
				my $cs      = $start;
				my $fstart  = $cur_hit->[2] > $cur_hit->[3] ? scalar(@frames) - 1 : 0;
				my $fdiff   = $cur_hit->[2] > $cur_hit->[3] ? -1 : 1;
				for ( my $fI = $fstart; $fI >= 0 && $fI < @frames; $fI += $fdiff ) {
					my $frame = $frames[$fI];
					if ( $frame =~ /(\d+):-*(\d+)/ ) {
						my $fs = $2;
						$fs *= -1 if $frame =~ /:-/;
						$cs += $fs;
					} else {
						$dna_seq .= substr( $seq, $cs - 1, $frame * 3 );
						$cs += $frame * 3;
					}
				}
				my $strand = 1;    # forward or reverse strand?
				$strand *= -1 if ( $cur_hit->[2] > $cur_hit->[3] );
				$new_seq = translate_frame(
											id                => $new_id,
											seq               => $dna_seq,
											frame             => $strand,
											marker            => $markerHit,
											reverse_translate => $self->{"dna"}
				);

				# bioperl uses * for stop codons but we want to give X to hmmer later
				$new_seq =~ s/\*/X/g;
			} elsif ( $type =~ /\.rna/ && $cur_hit->[2] > $cur_hit->[3] ) {

				# check if we need to reverse complement this one
				$new_seq =~ tr/ACGTacgt/TGCAtgca/;
				$new_seq = reverse($new_seq);
			}
			$markerHits{$markerHit} = "" unless defined( $markerHits{$markerHit} );
			$markerHits{$markerHit} .= ">".$new_id."\n".$new_seq."\n";
		}

		# ensure memory gets freed
		for ( my $i = 0; $i < $cur_hit_count; $i++ ) {
			$contig_hits{$id}->[$i] = undef;
		}
		$contig_hits{$id} = undef;
	}
	debug "$type Got ".scalar( keys %markerHits )." markers with hits\n";

	#write the read+ref_seqs for each markers in the list
	foreach my $marker ( keys %markerHits ) {

		#writing the hits to the candidate file
		my $candidate_file = Phylosift::Utilities::get_candidate_file(
																	   self   => $self,
																	   marker => $marker,
																	   type   => $type,
																	   chunk  => $chunk
		);
		my $FILE_OUT = ps_open( ">$candidate_file".".".$process_id );
		print $FILE_OUT $markerHits{$marker};
		close($FILE_OUT);
		if ( $self->{"dna"} && $type !~ /\.rna/ ) {
			$candidate_file = Phylosift::Utilities::get_candidate_file(
																		self   => $self,
																		marker => $marker,
																		type   => $type,
																		dna    => 1,
																		chunk  => $chunk
			);
			my $FILE_OUT = ps_open( ">$candidate_file".".".$process_id );
			print $FILE_OUT $markerNuc{$marker} if defined( $markerNuc{$marker} );
			close($FILE_OUT);
		}

		# force memory free
		$markerHits{$marker} = undef;
		$markerNuc{$marker}  = undef;
	}
	%markerNuc   = ();
	%contig_hits = ();
	%markerHits  = ();
}

=head2 prep_and_clean
=over

=item *

Checks if the directories needed for the blast run and parsing exist
Removes previous blast runs data if they are still in the directories
Generates the blastable database using the marker representatives

=back

=cut

sub prep_and_clean {
	my %args = @_;
	my $self = $args{self} || miss("self");
	`rm -rf "$self->{"blastDir"}"` if $Phylosift::Settings::force;

	#create a directory for the Reads file being processed.
	`mkdir "$Phylosift::Settings::file_dir"` unless ( -e $Phylosift::Settings::file_dir );
	`mkdir "$self->{"blastDir"}"` unless ( -e $self->{"blastDir"} );
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
