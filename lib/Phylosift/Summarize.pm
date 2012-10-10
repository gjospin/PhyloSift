package Phylosift::Summarize;
use warnings;
use strict;
use FindBin;
use Phylosift::Settings;
use Phylosift::Utilities qw(:all);
use Carp;
use Bio::Phylo;
use Bio::Phylo::Forest::Tree;
use IO::File;
use JSON;

set_default_values();

=head1 NAME

Phylosift::Summarize - Summarize placed reads using the NCBI taxonomy

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Reconciles gene-tree specific read placements with the NCBI taxonomy

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 summarize

=cut

my %nameidmap;
my %idnamemap;
my %ncbi_summary;    #used for krona output
my %deleted;
=head2 read_ncbi_taxon_name_map

# read the NCBI taxon names
# stash them in hashes called nameidmap and idnamemap to go back & forth from tax ids to names

=cut

sub read_ncbi_taxon_name_map {
	return ( \%nameidmap, \%idnamemap ) if %nameidmap;
	my $ncbidir = $Phylosift::Settings::ncbi_dir;
	my $TAXIDS  = ps_open("$ncbidir/names.dmp");
	while ( my $line = <$TAXIDS> ) {
		chomp $line;
		if (    ( $line =~ /scientific name/ )
			 || ( $line =~ /synonym/ )
			 || ( $line =~ /misspelling/ ) )
		{
			my @vals = split( /\s+\|\s+/, $line );
			$nameidmap{ homogenize_name_ala_dongying( name => $vals[1] ) } = $vals[0];
			$idnamemap{ $vals[0] } = homogenize_name_ala_dongying( name => $vals[1] )
			  if ( $line =~ /scientific name/ );
		}
	}
	return ( \%nameidmap, \%idnamemap );
}

# now read the NCBI taxonomy structure
# puts the results in a hash called "parent"
my %parent;

sub read_ncbi_taxonomy_structure {
	return \%parent if %parent;
	debug "Reading NCBI taxonomy\n";
	my $ncbidir      = $Phylosift::Settings::ncbi_dir;
	my $TAXSTRUCTURE = ps_open("$ncbidir/nodes.dmp");
	while ( my $line = <$TAXSTRUCTURE> ) {
		chomp $line;
		my @vals = split( /\s+\|\s+/, $line );
		$parent{ $vals[0] } = [ $vals[1], $vals[2] ];
	}
	return \%parent;
}

=head2 read_coverage
Reads a coverage file
Input: file - a file name
=cut

sub read_coverage {
	my %args = @_;
	my $file = $args{file} || miss("file");
	my %coverage;
	my $COVERAGE = ps_open($file);
	while ( my $line = <$COVERAGE> ) {
		chomp $line;
		my @data = split( /\t/, $line );
		$data[0] =~ s/[\(\)\[\]\+\=\<\>\?]//g;
		$data[0] =~ s/[\-\s]/_/g;
		$coverage{ $data[0] } = $data[1];
	}
	return \%coverage;
}

sub read_taxonmap {
	my %args     = @_;
	my $file     = $args{file} || miss("file");
	my $TAXONMAP = ps_open($file);
	my %markerncbimap;
	while ( my $line = <$TAXONMAP> ) {
		chomp($line);
		my ( $markerbranch, $ncbiname ) = split( /\t/, $line );
		$markerncbimap{$markerbranch} = []
		  unless defined( $markerncbimap{$markerbranch} );
		push( @{ $markerncbimap{$markerbranch} }, $ncbiname );
	}
	return \%markerncbimap;
}

=head2 set_default_values

set_default_values for all the parameters in this module

=cut

sub set_default_values {
	my %args = @_;
	my $self = $args{self};
	Phylosift::Settings::set_default(    parameter => \$Phylosift::Settings::krona_threshold,
									  value     => 0.01 );
}

=head2 summarize
Reads the .place files containing Pplacer read placements and maps them onto the
NCBI taxonomy
=cut

# keep a hash counting up all the read placements
# make this File-scope so anonymous functions below can see it
my %ncbireads;

sub summarize {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $chunk   = $args{chunk} || miss("chunk");
	my $markRef = $args{marker_reference} || miss("marker_reference");    # list of the markers we're using
	                                  #set_default_values(self=>$self);
	read_ncbi_taxon_name_map();
	read_ncbi_taxonomy_structure();
	my $markerdir = $Phylosift::Settings::marker_dir;

	# try to read a contig coverage file if it exists
	my %coverage;
	if ( defined $Phylosift::Settings::coverage ) {
		my $covref = read_coverage( file => $Phylosift::Settings::coverage );
		%coverage = %$covref;
	}

	# read all of the .place files for markers
	# map them onto the ncbi taxonomy
	# this is a hash structured as {sequenceID}{taxonID}=probabilitySum
	my %placements;
	my %unclassifiable;      # {sequenceID}=mass
	my %sequence_markers;    # {sequenceID}=[markers hit]
	my %weights;             #{sequenceID}{taxonID} = weight
	unshift( @{$markRef}, "concat" ) if $Phylosift::Settings::updated;
	for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {

		foreach my $marker ( @{$markRef} ) {

			# don't bother with this one if there's no read placements
			my $sub_mark;
			$sub_mark = "*" if $dna;
			my $place_base = $self->{"treeDir"} . "/"
			  . Phylosift::Utilities::get_read_placement_file( marker => $marker, dna => $dna, sub_marker => $sub_mark, chunk => $chunk );
			my @place_files;
			@place_files = glob($place_base) if $dna;    # need to glob on all submarkers if in DNA
			push( @place_files, $place_base ) if !$dna && -e $place_base;
			foreach my $placeFile (@place_files) {
				my $PP_COVFILE = ps_open(">".Phylosift::Utilities::get_read_placement_file(marker => $marker, chunk  => $chunk). ".cov") if ( defined $Phylosift::Settings::coverage );
				my $sub;
				$sub = $1 if $placeFile =~ /\.sub(\d+)\./;

				# first read the taxonomy mapping
				my $markermapfile = Phylosift::Utilities::get_marker_taxon_map(self => $self, marker => $marker, dna => $dna, sub_marker => $sub);
				next unless -e $markermapfile;    # can't summarize if there ain't no mappin'!
				my $markerncbimap = read_taxonmap( file => $markermapfile );

				# then read & map the placement
				my $JPLACEFILE = ps_open($placeFile);
				my @treedata   = <$JPLACEFILE>;
				close $JPLACEFILE;
				my $json_data = decode_json( join( "", @treedata ) );

				# for each placement record
				for ( my $i = 0 ; $i < @{ $json_data->{placements} } ; $i++ ) {
					my $place = $json_data->{placements}->[$i];

					# for each placement edge in the placement record
					for ( my $j = 0 ; $j < @{ $place->{p} } ; $j++ ) {
						my $edge      = $place->{p}->[$j]->[0];
						my $edge_mass = $place->{p}->[$j]->[2];
						if ( !defined( $markerncbimap->{$edge} ) ) {

							# mark these reads as unclassifiable
							for ( my $k = 0 ; $k < @{ $place->{nm} } ; $k++ ) {
								my $qname   = $place->{nm}->[$k]->[0];
								my $qweight = $place->{nm}->[$k]->[1];
								$unclassifiable{$qname} = 0 unless defined( $unclassifiable{$qname} );
								$unclassifiable{$qname} += $edge_mass * $qweight;
							}
							next;
						}
						my $mapcount = scalar( @{ $markerncbimap->{$edge} } );

						# for each taxon to which the current phylogeny edge could map
						foreach my $taxon ( @{ $markerncbimap->{$edge} } ) {
							my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
							
							# for each query seq in the current placement record
							for ( my $k = 0 ; $k < @{ $place->{nm} } ; $k++ ) {
								my $qname   = $place->{nm}->[$k]->[0];
								my $qweight = $place->{nm}->[$k]->[1];
								$sequence_markers{$qname}{$marker} = 1;
								$placements{$qname} = () unless defined( $placements{$qname} );
								$placements{$qname}{$taxon_id} = 0 unless defined( $placements{$qname}{$taxon_id} );
								$placements{$qname}{$taxon_id} += $edge_mass * $qweight / $mapcount;
								$ncbireads{$taxon} = 0 unless defined $ncbireads{$taxon};

								# split the p.p. across the possible edge mappings
								$ncbireads{$taxon} += $edge_mass * $qweight / $mapcount;
							}
						}
					}

				}
			}
		}

	}

	# also write out the taxon assignments for sequences
	write_sequence_taxa_summary( self=>$self, placements => \%placements, unclassifiable => \%unclassifiable, sequence_markers => \%sequence_markers, chunk => $chunk );
	Phylosift::Summarize::rename_sequences( self => $self, chunk => $chunk );
	$self->{"read_names"} = ();    #clearing the hash from memory
	merge_sequence_taxa( self => $self, chunk => $chunk );
}

sub write_sequence_taxa_summary {
	my %args  = @_;
	my $self  = $args{self} || miss("PS Object");
	my $placements = $args{placements} || miss("placements");
	my $unclassifiable = $args{unclassifiable} || miss("unclassifiable");
	my $sequence_markers = $args{sequence_markers} || miss("sequence_markers");
	my $chunk = $args{chunk};
	my $chunky = defined($chunk) ? ".$chunk" : "";
	
	my $SEQUENCETAXA = ps_open( ">$Phylosift::Settings::file_dir/sequence_taxa$chunky.txt" );
	my $SEQUENCESUMMARY = ps_open( ">$Phylosift::Settings::file_dir/sequence_taxa_summary$chunky.txt" );
	foreach my $qname ( keys(%$placements) ) {

		# sum up all placements for this sequence, use to normalize
		my $placecount = 0;
		foreach my $taxon_id ( keys( %{ $placements->{$qname} } ) ) {
			$placecount += $placements->{$qname}{$taxon_id};
		}
		$placecount += $unclassifiable->{$qname} if defined( $unclassifiable->{$qname} );

		# determine the different unique names used for this molecule (e.g. /1 and /2 for paired reads)
		my %unique_names;

		# normalize to probability distribution
		foreach my $taxon_id ( sort { $placements->{$qname}{$b} <=> $placements->{$qname}{$a} } keys %{ $placements->{$qname} } )
		{
			#skip if taxon_id is empty/undefined
			next unless defined($taxon_id) || $taxon_id eq "";
			my $blah = $taxon_id || next;
			$placements->{$qname}{$taxon_id} /= $placecount;
			my ( $taxon_name, $taxon_level, $tid ) =  get_taxon_info( taxon => $taxon_id );
			$taxon_level = "Unknown" unless defined($taxon_level);
			$taxon_name  = "Unknown" unless defined($taxon_name);
			$self->{"read_names"}{$qname} = [$qname] unless defined( $self->{"read_names"}{$qname} );
			if ( exists $self->{"read_names"}{$qname} ) {

				#				$placements{$qname}{$taxon_id} /=
				#				  @{ $self->{"read_names"}{$qname} };
				foreach my $name_ref ( @{ $self->{"read_names"}{$qname} } ) {
					$unique_names{$name_ref} = 1;
				}
				foreach my $name_ref ( keys(%unique_names) ) {
					print $SEQUENCETAXA "$name_ref\t$taxon_id\t$taxon_level\t$taxon_name\t"
					  . $placements->{$qname}{$taxon_id} . "\t" . join( "\t", keys( %{ $sequence_markers->{$qname} } ) ). "\n";
				}
			}
		}
		foreach my $name_ref ( keys(%unique_names) ) {
			if ( defined( $unclassifiable->{$qname} ) ) {
				print $SEQUENCETAXA "$name_ref\tUnknown\tUnknown\tUnclassifiable\t" . ( $unclassifiable->{$qname} / $placecount ) . "\t"
				  . $unclassifiable->{$qname} . "\t" . join( "\t", keys( %{ $sequence_markers->{$qname} } ) ) . "\n";
			}
		}

		my $readsummary = sum_taxon_levels( placements => $placements->{$qname} );
		foreach my $taxon_id ( sort { $readsummary->{$b} <=> $readsummary->{$a} } keys %{$readsummary} ) {
			my $blah = $taxon_id || next;
			my ( $taxon_name, $taxon_level, $tid ) = get_taxon_info( taxon => $taxon_id );
			$taxon_level = "Unknown" unless defined($taxon_level);
			$taxon_name  = "Unknown" unless defined($taxon_name);
			if ( exists $self->{"read_names"}{$qname} ) {
				foreach my $name_ref ( keys(%unique_names) ) {
					print $SEQUENCESUMMARY "$name_ref\t$taxon_id\t$taxon_level\t$taxon_name\t" . $readsummary->{$taxon_id} . "\t"
					  . join( "\t", keys( %{ $sequence_markers->{$qname} } ) ) . "\n";
				}
			}

		}
	}
	close($SEQUENCESUMMARY);
	close($SEQUENCETAXA);
}

=head2 merge_sequence_taxa

	Reads all the sequence_taxa files and outputs 1 single taxasummary and other statistics files for the whole sample.

=cut

sub merge_sequence_taxa {
	my %args       = @_;
	my $self       = $args{self} || miss("PS Object");
	my $chunk      = $args{chunk} || miss("chunk");
	my $taxa_seed  = $Phylosift::Settings::file_dir . "/sequence_taxa.*.txt";
	my @taxa_files = glob("$taxa_seed");
	my %placements = ();
	my %placement_markers = ();
	my %all_summary;
	my %concat_summary;
	my %unclassifiable;    # {sequenceID}=mass
	my $unclass_total          = 0;
	my %marker_number_hits     = ();
	my %unclassifiable_markers = ();

	foreach my $taxa_file (@taxa_files) {

		#read all the sequence information
		my $TAXAIN = ps_open($taxa_file);
		while (<$TAXAIN>) {
			chomp($_);
			my @line         = split( /\t/, $_ );
			my $read_id      = $line[0];
			my $taxon_id     = $line[1];
			my $rank         = $line[2];
			my $species_name = $line[3];
			my $prob         = $line[4];
			my $marker_name  = $line[5];
			if ( $taxon_id eq "Unknown" ) {
				$unclassifiable{$read_id} = $prob;
				for ( my $i = 5 ; $i < @line ; $i++ ) {
					$unclassifiable_markers{$read_id}{ $line[$i] } = 1;
				}
			} else {
				$placements{$read_id}{$taxon_id} = $prob;
				for ( my $i = 5 ; $i < @line ; $i++ ) {
					$placement_markers{$read_id}{ $line[$i] } = 1;
				}
			}
		}

		# make a summary of total reads at each taxonomic level
		# this gets used later in krona output
		foreach my $qname ( keys(%placements) ) {
			my $readsummary = sum_taxon_levels( placements => $placements{$qname} );
			foreach my $taxon_id ( keys(%$readsummary) ) {
				$all_summary{$taxon_id} = 0 unless defined( $all_summary{$taxon_id} );
				$all_summary{$taxon_id} += $readsummary->{$taxon_id};
				next unless defined( $placement_markers{$qname}{"concat"} );
				$concat_summary{$taxon_id} = 0 unless defined( $concat_summary{$taxon_id} );
				$concat_summary{$taxon_id} += $readsummary->{$taxon_id};
			}
		}

		#gathering hit numbers per markers
		foreach my $rid ( keys %unclassifiable_markers ) {
			foreach my $mname ( keys %{ $unclassifiable_markers{$rid} } ) {
				if ( exists $marker_number_hits{$mname} ) {
					$marker_number_hits{$mname}++;
				} else {
					$marker_number_hits{$mname} = 1;
				}
			}
		}
		foreach my $rid ( keys %placement_markers ) {
			foreach my $mname ( keys %{ $placement_markers{$rid} } ) {
				if ( exists $marker_number_hits{$mname} ) {
					$marker_number_hits{$mname}++;
				} else {
					$marker_number_hits{$mname} = 1;
				}
			}
		}

		# write unclassifiable

		foreach my $u ( values(%unclassifiable) ) {
			$unclass_total += $u;
		}

		close($TAXAIN);
		%placements             = ();
		%placement_markers      = ();
		%unclassifiable         = ();
		%unclassifiable_markers = ();
	}

	# write a single combined taxa summary
	`cat $taxa_seed > $Phylosift::Settings::file_dir/sequence_taxa.txt`;
	my $taxa_summary_seed  = $Phylosift::Settings::file_dir . "/sequence_taxa_summary.*.txt";
	`cat $taxa_summary_seed > $Phylosift::Settings::file_dir/sequence_taxa_summary.txt`;


	my $MARKER_HITS = ps_open( ">" . $Phylosift::Settings::file_dir . "/marker_summary.txt" );
	print $MARKER_HITS "Marker_name\tNumber of hits\n";
	foreach my $mname ( sort { $marker_number_hits{$b} <=> $marker_number_hits{$a} } keys %marker_number_hits ) {
		print $MARKER_HITS "$mname\t$marker_number_hits{$mname}\n";
	}
	close($MARKER_HITS);
	my $TAXAOUT = ps_open( ">" . $Phylosift::Settings::file_dir . "/taxasummary.txt" );

	# sort rest of taxa by descending abundance order
	print $TAXAOUT "Unclassifiable\tUnknown\tUnknown\t$unclass_total\n";
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		$taxon_level = "Unknown" unless defined($taxon_level);
		$taxon_name  = "Unknown" unless defined($taxon_name);
		print $TAXAOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
	}
	close($TAXAOUT);

	# sample from multinomial to get confidence limits
	# get total read count
	my $totalreads = 0;
	foreach my $val ( values(%ncbireads) ) {
		$totalreads += $val;
	}
	debug "Total reads placed is " . scalar( keys(%placements) ) . "\n";
	debug "Total classifiable probability mass is $totalreads\n";

	# write the taxa with 90% highest posterior density, assuming each read is an independent observation
	my $taxasum    = 0;
	my $TAXAHPDOUT = ps_open( ">" . $Phylosift::Settings::file_dir . "/taxa_90pct_HPD.txt" );
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		$taxasum += $ncbireads{$taxon};
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		print $TAXAHPDOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
		last if $taxasum >= $totalreads * 0.9;
	}
	close($TAXAHPDOUT);
	Phylosift::Utilities::start_timer( name => "runKrona" );

	#Need to move this to the merge summary function
	unless ($Phylosift::Settings::simple) {

		# skip this if only a simple summary is desired (it's slow)
		debug "Generating krona\n";
		my $html_report = $Phylosift::Settings::file_dir."/".$self->{"fileName"}.".html";
		#$self->{HTML} = ps_open(">$html_report") unless defined $self->{HTML};
		if(scalar( keys(%concat_summary) ) > 0){
			if(!defined $self->{HTML} ){
				$self->{HTML} = Phylosift::HTMLReport::begin_report( self=>$self, file => $html_report );
			}
			Phylosift::HTMLReport::add_krona(self=>$self, OUTPUT=>$self->{HTML}, summary => \%concat_summary );
		}else{
			my $FH = ps_open(">$html_report");
			print $FH "No results Found\n";
			close($FH);
		}
		%ncbi_summary = ();
		%ncbi_summary = %all_summary;
		my $html_all_report = $Phylosift::Settings::file_dir."/".$self->{"fileName"}."_allmarkers.html";
		if(scalar( keys(%ncbi_summary) ) > 0){
			if(!defined $self->{HTMLall} ){
				$self->{HTMLall} = Phylosift::HTMLReport::begin_report( self=>$self, file => $html_all_report );
			}
			Phylosift::HTMLReport::add_krona(self=>$self, OUTPUT=>$self->{HTMLall}, summary => \%ncbi_summary );
		}else{
			my $FH = ps_open(">$html_all_report");
			print $FH "No results Found\n";
			close($FH);
		}
	}
	Phylosift::Utilities::end_timer( name => "runKrona" );

}

sub print_run_info {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $OUTPUT  = $args{OUTPUT} || miss("OUTPUT");
	my $newline = $args{newline} || "\n";
	print $OUTPUT "Program run as:$newline" . join( " ", "phylosift", @{ $self->{"ARGV"} } ) . "$newline";
	print $OUTPUT "Marker database version:\n" . join("$newline", 
		Phylosift::Utilities::get_marker_version( path => $Phylosift::Settings::marker_dir ) ). "$newline";
}


#
# non-functional until dependency on Math::Random can be eliminated
#
sub write_confidence_intervals {
	my %args         = @_;
	my $self         = $args{self} || miss("self");
	my $ncbireadsref = $args{ncbi_reads_reference} || miss("ncbi_reads_reference");
	my $totalreads = $args{total_reads} || miss("total_reads");
	my %ncbireads  = %$ncbireadsref;

	# normalize to a sampling distribution
	foreach my $key ( keys(%ncbireads) ) {
		$ncbireads{$key} /= $totalreads + 1;
	}
	my $normsum  = 0;
	my @valarray = values(%ncbireads);
	foreach my $val (@valarray) {
		$normsum += $val;
	}
	my $sample_count = 100;
	my %samples;
	for ( my $sI = 0 ; $sI < $sample_count ; $sI++ ) {

		#        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
		my @sample;
		my $kI = 0;
		foreach my $key ( keys(%ncbireads) ) {
			push( @{ $samples{$key} }, $sample[ $kI++ ] );
		}
	}
	my $TAXA_CONF = ps_open( ">" . $Phylosift::Settings::file_dir . "/taxaconfidence.txt" );
	foreach my $key ( keys(%samples) ) {
		my @svals = @{ $samples{$key} };
		my @sorted = sort { $a <=> $b } @svals;
		my ( $taxon_name, $taxon_level, $taxon_id ) = getTaxonInfo($key);
		print $TAXA_CONF join( "\t",
							   $taxon_id,
							   $taxon_level,
							   $taxon_name,
							   $sorted[0],
							   $sorted[ int( $sample_count * 0.1 ) ],
							   $sorted[ int( $sample_count * 0.25 ) ],
							   $sorted[ int( $sample_count * 0.5 ) ],
							   $sorted[ int( $sample_count * 0.75 ) ],
							   $sorted[ int( $sample_count * 0.9 ) ],
							   $sorted[ $sample_count - 1 ] ),
		  "\n";
	}
}

sub sum_taxon_levels {
	my %args       = @_;
	my $placements = $args{placements} || miss("placements");
	my %summarized = ();
	foreach my $taxon_id ( keys %$placements ) {
		my $cur_tid = $taxon_id;
		while ( defined($cur_tid) ) {
			$summarized{$cur_tid} = 0 unless defined( $summarized{$cur_tid} );
			$summarized{$cur_tid} += $placements->{$taxon_id};
			last if defined( $parent{$cur_tid}[0] ) && $parent{$cur_tid}[0] == $cur_tid;
			$cur_tid = $parent{$cur_tid}[0];
		}
	}
	return \%summarized;
}

sub get_taxon_info {
	my %args = @_;
	my $in = $args{taxon} || miss("taxon");
	read_ncbi_taxon_name_map;
	read_ncbi_taxonomy_structure();
	if ( $in =~ /^\d+$/ ) {

		#it's an ncbi taxon id.  look up its name and level.
		my $merged = Phylosift::Summarize::read_merged_nodes();
		$in = $merged->{$in} if defined( $merged->{$in} );
		my $deleted = Phylosift::Summarize::read_deleted_nodes();
		return ($in, "", "") if defined($deleted->{$in});
		my $name  = $idnamemap{$in};
		my $level = $parent{$in}->[1];
		return ( $name, $level, $in );
	} elsif ( $in =~ /\w+/ ) {

		# old style map, need to go from NCBI name back to ID
		my $name = $in;
		$name =~ s/_/ /g;    # map uses spaces instead of underscores
		my ( $id, $qname ) = donying_find_name_in_taxa_db( name => $name );
		my $level = $parent{$id}->[1];
		return ( $in, $level, $id );
	}
	return ( $in, "", "" );
}

sub tree_name {
	my %args = @_;
	my $inName = $args{name} || return;
	$inName =~ s/\s+/_/g;
	$inName =~ s/'//g;
	$inName =~ s/[\(\)]//g;
	$inName =~ s/-/_/g;
	$inName =~ s/\//_/g;
	$inName =~ s/#/_/g;
	$inName =~ s/\:/_/g;
	return $inName;
}

=head2 homogenize_name_ala_dongying

=cut

sub homogenize_name_ala_dongying {
	my %args = @_;
	my $inName = $args{name} || miss("name");
	return "" unless defined($inName);
	$inName =~ s/^\s+//;
	$inName =~ s/\s+$//;
	$inName =~ s/\s+/ /g;
	$inName =~ s/,//g;
	$inName =~ s/[)(]//g;
	$inName =~ s/-/ /g;
	$inName = uc $inName;
	return $inName;
}

=head2 donying_find_name_in_taxa_db
    
=cut

sub donying_find_name_in_taxa_db {
	my %args = @_;
	my $name = $args{name} || miss("name");
	return "" unless defined($name);
	$name =~ s/^\s+//;
	my @t          = split( /\s+/, $name );
	my $input_name = join( " ",    @t );
	my $q_name     = $input_name;
	my $id         = "ERROR";
	while ( @t >= 1 ) {
		$q_name = join( " ", @t );
		$q_name = uc $q_name;
		if ( defined( $nameidmap{$q_name} ) ) {
			$id = $nameidmap{$q_name};
			last;
		}
		pop @t;
	}
	return ( $id, $q_name );
}

=head2 read_merged_nodes

Read in obsolete NCBI taxon IDs that have been merged into new taxon IDs.
Returns a hash mapping old to new taxon ID

=cut

my %merged;

sub read_merged_nodes {
	return \%merged if %merged;
	debug "Reading merged ncbi nodes\n";
	my $MERGED = ps_open("$Phylosift::Settings::ncbi_dir/merged.dmp");
	while ( my $line = <$MERGED> ) {
		chomp $line;
		my @vals = split( /\s+\|\s*/, $line );
		$merged{ $vals[0] } = $vals[1];
	}
	debug "Done reading merged\n";
	return \%merged;
}

=head2 rename_sequences

Looks through all the output files in the directories to change the unique numbers back to the original sequence IDs
Only need a PS object as input.
=cut

sub rename_sequences {
	my %args  = @_;
	my $self  = $args{self} || miss("PS object");
	my $chunk = $args{chunk} || miss("Chunk number");
	if(!-e $self->{"blastDir"} . "/lookup_ID.$chunk.tbl"){
		warn "Lookup_ID.$chunk.tbl was not found\n";
		return;
	}
	#read the read mapping file
	my %name_mapping = ();
	my $FH_MAP       = ps_open( $self->{"blastDir"} . "/lookup_ID.$chunk.tbl" );
	while (<$FH_MAP>) {
		chomp($_);
		my @line = split( / /, $_ );
		$name_mapping{ $line[1] } = $line[0];
	}
	close($FH_MAP);

	#open all the
	my @array_to_rename = ();
	push( @array_to_rename, glob( $self->{"blastDir"} . "/*.ffn.$chunk*" ) );
	push( @array_to_rename, glob( $self->{"blastDir"} . "/*.aa.$chunk*" ) );
	push( @array_to_rename, glob( $self->{"alignDir"} . "/*.updated.$chunk.fasta" ) );
	push( @array_to_rename, glob( $self->{"alignDir"} . "/*.newCandidate.aa.$chunk" ) );
	push( @array_to_rename, glob( $self->{"alignDir"} . "/*.$chunk.unmasked" ) );
	push( @array_to_rename, glob( $self->{"alignDir"} . "/*.codon.updated.$chunk.fasta" ) );
	push( @array_to_rename, glob( $self->{"treeDir"} . "/*.updated.$chunk.jplace" ) );
	push( @array_to_rename, glob( $Phylosift::Settings::file_dir . "/*.jplace" ) );
	push( @array_to_rename, glob( $Phylosift::Settings::file_dir . "/sequence_taxa*.$chunk.txt" ) );
	debug "FILE DIR : ".$Phylosift::Settings::file_dir ."\n";
	foreach my $file (@array_to_rename) {
		my $FH  = ps_open($file);
		my $TMP = ps_open(">$file.tmp");
		if ( $file =~ m/\.jplace/ ) {
			my @treedata = <$FH>;	
			my $json_data = decode_json( join( "", @treedata ) );
			# parse the tree
			for ( my $i = 0 ; $i < @{ $json_data->{placements} } ; $i++ ) {
				my $placement = $json_data->{placements}->[$i];
				for ( my $j = 0 ; $j < @{ $placement->{nm} } ; $j++ ) {
					$placement->{nm}->[$j]->[0] =~ s/^(\d+)\./$name_mapping{$1}\./;
				}
			}

			# write the renamed jplace
			print $TMP encode_json($json_data);	
		}else{
			while (<$FH>) {
				$_ =~ s/^>(\d+)\./>$name_mapping{$1}\./g; #fasta files
				$_ =~ s/^(\d+)\./$name_mapping{$1}\./g; #summary files
				print $TMP $_;
			}
		}
		close($FH);
		close($TMP);
		`mv "$file.tmp" "$file"`;
	}
}

sub read_deleted_nodes {
	return \%deleted if %deleted;
	debug "Reading deleted ncbi nodes\n";
	my $DELETED = ps_open("$Phylosift::Settings::ncbi_dir/delnodes.dmp");
	while ( my $line = <$DELETED> ) {
		chomp $line;
		my @vals = split( /\s+\|\s*/, $line );
		$deleted{ $vals[0] } = 1;
	}
	debug "Done reading deleted\n";
	return \%deleted;
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

    perldoc Phylosift::Summarize


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

1;    # End of Phylosift::Summarize.pm
