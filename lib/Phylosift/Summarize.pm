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
use File::Basename;
our $VERSION = "v1.0.1";

set_default_values();

=head1 NAME

Phylosift::Summarize - Summarize placed reads using the NCBI taxonomy

=head1 VERSION

Version 0.01

=cut

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
	debug "Reading NCBI taxonomy at $Phylosift::Settings::ncbi_dir\n";
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
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::krona_threshold,
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

	Phylosift::Utilities::start_step( self => $self, chunk => $chunk, step => "Summarize" );
	my $completed_chunk = Phylosift::Utilities::has_step_completed( self => $self, chunk => $chunk, step => "Summarize", force => $Phylosift::Settings::force );
	croak("Previous step for chunk $chunk has did not complete. Aborting\n")
	  unless Phylosift::Utilities::has_step_completed( self => $self, chunk => $chunk, step => "Place" );
	unless ($completed_chunk) {

		#read_ncbi_taxon_name_map();
		#read_ncbi_taxonomy_structure();
		my $markerdir = $Phylosift::Settings::marker_dir;

		# try to read a contig coverage file if it exists
		my %coverage;
		if ( defined $Phylosift::Settings::coverage ) {
			my $covref = read_coverage( file => $Phylosift::Settings::coverage );
			%coverage = %$covref;
		}
		debug "\n\n\n******STARTING SUMMARY \n\n\n";

		# read all of the .place files for markers
		# map them onto the ncbi taxonomy
		# this is a hash structured as {sequenceID}{taxonID}=probabilitySum
		my %placements;
		my %unclassifiable;      # {sequenceID}=mass
		my %sequence_markers;    # {sequenceID}=[markers hit]
		my %weights;             #{sequenceID}{taxonID} = weight
		                         #unshift( @{$markRef}, "concat" ) if $Phylosift::Settings::updated;
		foreach my $marker ( @{$markRef} ) {

			my %seen_queries;    # tracks whether a query has been seen already
			for ( my $dna = 1; $dna >= 0; $dna-- ) {
				my $sub = $dna ? 1 : undef;
				my $place_file =
				  $self->{"treeDir"}."/"
				  .Phylosift::Utilities::get_read_placement_file( self => $self, marker => $marker, dna => $dna, sub_marker => $sub, chunk => $chunk );
				next unless -e $place_file;    # don't bother with this one if there's no read placements
				my $PP_COVFILE = ps_open( ">".Phylosift::Utilities::get_read_placement_file( marker => $marker, chunk => $chunk ).".cov" )
				  if ( defined $Phylosift::Settings::coverage );

				# first read the taxonomy mapping
				my $markermapfile = Phylosift::Utilities::get_marker_taxon_map( self => $self, marker => $marker, dna => $dna );
				debug "Marker $marker missing taxon map $markermapfile\n" unless -e $markermapfile;
				next unless -e $markermapfile;    # can't summarize if there ain't no mappin'!
				my $markerncbimap = read_taxonmap( file => $markermapfile );

				# then read & map the placement
				my $JPLACEFILE = ps_open($place_file);
				my @treedata   = <$JPLACEFILE>;
				close $JPLACEFILE;
				my $json_data = decode_json( join( "", @treedata ) );

				# for each placement record
				for ( my $i = 0; $i < @{ $json_data->{placements} }; $i++ ) {
					my $place = $json_data->{placements}->[$i];

					# for each query in the placement record
					for ( my $k = 0; $k < @{ $place->{nm} }; $k++ ) {
						my $qname   = $place->{nm}->[$k]->[0];
						my $qweight = $place->{nm}->[$k]->[1];

						# skip this query if it was already processed (e.g. in codons)
						next if defined( $seen_queries{$qname} );
						$seen_queries{$qname} = 1;

						# for each placement edge in the placement record
						for ( my $j = 0; $j < @{ $place->{p} }; $j++ ) {
							my $edge      = $place->{p}->[$j]->[1];
							my $edge_mass = $place->{p}->[$j]->[2];
							if ( !defined( $markerncbimap->{$edge} ) ) {

								# mark these reads as unclassifiable
								$unclassifiable{$qname} = 0 unless defined( $unclassifiable{$qname} );
								$unclassifiable{$qname} += $edge_mass * $qweight;
								next;
							}

							# for each taxon to which the current phylogeny edge could map
							# add the placement probability mass for the sequence to that taxonomic group
							my $mapcount = scalar( @{ $markerncbimap->{$edge} } );
							foreach my $taxon ( @{ $markerncbimap->{$edge} } ) {
								my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );

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

		# also write out the taxon assignments for sequences
		write_sequence_taxa_summary(
									 self             => $self,
									 placements       => \%placements,
									 unclassifiable   => \%unclassifiable,
									 sequence_markers => \%sequence_markers,
									 chunk            => $chunk
		);
		Phylosift::Summarize::rename_sequences( self => $self, chunk => $chunk );
		$self->{"read_names"} = ();    #clearing the hash from memory
		merge_sequence_taxa( self => $self, chunk => $chunk );
	}
	Phylosift::Utilities::end_step( self => $self, chunk => $chunk, step => "Summarize" );
	Phylosift::Utilities::write_step_completion_to_run_info( self => $self, chunk => $chunk, step => "Summarize" ) unless $completed_chunk;
}

sub write_sequence_taxa_summary {
	my %args             = @_;
	my $self             = $args{self} || miss("PS Object");
	my $placements       = $args{placements} || miss("placements");
	my $unclassifiable   = $args{unclassifiable} || miss("unclassifiable");
	my $sequence_markers = $args{sequence_markers} || miss("sequence_markers");
	my $chunk            = $args{chunk};
	my $chunky           = defined($chunk) ? ".$chunk" : "";
	debug "Writing sequences\n";
	my $SEQUENCETAXA = ps_open(">$Phylosift::Settings::file_dir/sequence_taxa$chunky.txt");
	print $SEQUENCETAXA "##Sequence_ID\tHit_Coordinates\tNCBI_Taxon_ID\tTaxon_Rank\tTaxon_Name\tProbability_Mass\tMarkers_Hit\n";
	my $SEQUENCESUMMARY = ps_open(">$Phylosift::Settings::file_dir/sequence_taxa_summary$chunky.txt");
	print $SEQUENCESUMMARY "#Sequence_ID\tHit_Coordinates\tNCBI_Taxon_ID\tTaxon_Rank\tTaxon_Name\tCumulative_Probability_Mass\tMarkers_Hit\n";

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
		foreach my $taxon_id ( sort { $placements->{$qname}{$b} <=> $placements->{$qname}{$a} } keys %{ $placements->{$qname} } ) {

			#skip if taxon_id is empty/undefined
			next unless defined($taxon_id) || $taxon_id eq "";
			my $blah = $taxon_id || next;
			$placements->{$qname}{$taxon_id} /= $placecount;
			my ( $taxon_name, $taxon_level, $tid ) = get_taxon_info( taxon => $taxon_id );
			$taxon_level                  = "Unknown" unless defined($taxon_level);
			$taxon_name                   = "Unknown" unless defined($taxon_name);
			$self->{"read_names"}{$qname} = [$qname]  unless defined( $self->{"read_names"}{$qname} );
			if ( exists $self->{"read_names"}{$qname} ) {

				#				$placements{$qname}{$taxon_id} /=
				#				  @{ $self->{"read_names"}{$qname} };
				foreach my $name_ref ( @{ $self->{"read_names"}{$qname} } ) {
					$unique_names{$name_ref} = 1;
				}
				foreach my $name_ref ( keys(%unique_names) ) {
					my %coords = ();
					my @coord_split = split( /\./, $name_ref );
					$name_ref = shift(@coord_split);    #removes the first element of the array
					for ( my $i = 0; $i < @coord_split; $i = $i + 3 ) {
						my $start = $coord_split[ $i + 1 ];
						my $end   = $coord_split[ $i + 2 ];
						if ( defined( $coords{ $coord_split[$i] } ) ) {
							$coords{ $coord_split[$i] } .= ".$start.$end";
						} else {
							$coords{ $coord_split[$i] } = "$start.$end";
						}
					}
					foreach my $mate ( keys %coords ) {
						print $SEQUENCETAXA "$name_ref/$mate\t$coords{$mate}\t$taxon_id\t$taxon_level\t$taxon_name\t"
						  .$placements->{$qname}{$taxon_id}."\t"
						  .join( "\t", keys( %{ $sequence_markers->{$qname} } ) )."\n";
					}
					%coords      = ();
					@coord_split = ();
				}
			}
		}
		foreach my $name_ref ( keys(%unique_names) ) {
			if ( defined( $unclassifiable->{$qname} ) ) {
				my %coords = ();
				my @coord_split = split( /\./, $name_ref );
				$name_ref = shift(@coord_split);    #removes the first element of the array
				for ( my $i = 0; $i < @coord_split; $i = $i + 3 ) {
					if ( defined( $coords{ $coord_split[$i] } ) ) {
						$coords{ $coord_split[$i] } .= ".$coord_split[$i+1].$coord_split[$i+2]";
					} else {
						$coords{ $coord_split[$i] } = "$coord_split[$i+1].$coord_split[$i+2]";
					}
				}
				foreach my $mate ( keys %coords ) {
					print $SEQUENCETAXA "$name_ref/$mate\t$coords{$mate}\tUnknown\tUnknown\tUnclassifiable\t".( $unclassifiable->{$qname} / $placecount )."\t"
					  ##.$unclassifiable->{$qname}."\t"
					  .join( "\t", keys( %{ $sequence_markers->{$qname} } ) )."\n";
				}
				%coords      = ();
				@coord_split = ();
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
					my %coords = ();
					my @coord_split = split( /\./, $name_ref );
					$name_ref = shift(@coord_split);    #removes the first element of the array
					for ( my $i = 0; $i < @coord_split; $i = $i + 3 ) {
						if ( defined( $coords{ $coord_split[$i] } ) ) {
							$coords{ $coord_split[$i] } .= ".$coord_split[$i+1].$coord_split[$i+2]";
						} else {
							$coords{ $coord_split[$i] } = "$coord_split[$i+1].$coord_split[$i+2]";
						}
					}
					foreach my $mate ( keys %coords ) {
						print $SEQUENCESUMMARY "$name_ref/$mate\t$coords{$mate}\t$taxon_id\t$taxon_level\t$taxon_name\t"
						  .$readsummary->{$taxon_id}."\t"
						  .join( "\t", keys( %{ $sequence_markers->{$qname} } ) )."\n";
					}

					%coords      = ();
					@coord_split = ();
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
	my %args              = @_;
	my $self              = $args{self} || miss("PS Object");
	my $chunk             = $args{chunk} || miss("chunk");
	my @taxa_files        = Phylosift::Utilities::get_summarize_output_sequence_taxa( self => $self );
	my %placements        = ();
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
			next if $_ =~ m/^#/;
			my @line         = split( /\t/, $_ );
			my $read_id      = $line[0];
			my $coord        = $line[1];
			my $taxon_id     = $line[2];
			my $rank         = $line[3];
			my $species_name = $line[4];
			my $prob         = $line[5];
			my $marker_name  = $line[6];
			if ( $taxon_id eq "Unknown" ) {
				$unclassifiable{$read_id} = $prob;
				for ( my $i = 6; $i < @line; $i++ ) {
					$unclassifiable_markers{$read_id}{ $line[$i] } = 1;
				}
			} else {
				$placements{$read_id}{$taxon_id} = $prob;
				for ( my $i = 6; $i < @line; $i++ ) {
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
	my @taxa_summary_files = Phylosift::Utilities::get_summarize_output_sequence_taxa_summary( self => $self );
	Phylosift::Utilities::concatenate_summary_files( self => $self, files => \@taxa_files, output_file => "$Phylosift::Settings::file_dir/sequence_taxa.txt" );
	Phylosift::Utilities::concatenate_summary_files(
													 self        => $self,
													 files       => \@taxa_summary_files,
													 output_file => "$Phylosift::Settings::file_dir/sequence_taxa_summary.txt"
	);

	my $MARKER_HITS = ps_open( ">".$Phylosift::Settings::file_dir."/marker_summary.txt" );
	print $MARKER_HITS "#Marker_name\tNumber of hits\n";
	foreach my $mname ( sort { $marker_number_hits{$b} <=> $marker_number_hits{$a} } keys %marker_number_hits ) {
		print $MARKER_HITS "$mname\t$marker_number_hits{$mname}\n";
	}
	close($MARKER_HITS);
	my $TAXAOUT = ps_open( ">".Phylosift::Utilities::get_taxasummary( self => $self ) );
	print $TAXAOUT "#Taxon_ID\tTaxon_Rank\tTaxon_Name\tProbability_Mass\n";

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
	debug "Total classifiable probability mass is $totalreads\n";

	# write the taxa with 90% highest posterior density, assuming each read is an independent observation
	my $taxasum = 0;
	my $TAXAHPDOUT = ps_open( ">".Phylosift::Utilities::get_taxa_90pct_HPD( self => $self ) );
	print $TAXAHPDOUT "#Taxon_ID\tTaxon_Rank\tTaxon_Name\tProbability_Mass\n";
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		$taxasum += $ncbireads{$taxon};
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		$taxon_level = "" unless defined $taxon_level;
		$taxon_name = "" unless defined $taxon_name;
		print $TAXAHPDOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
		last if $taxasum >= $totalreads * 0.9;
	}
	close($TAXAHPDOUT);
	Phylosift::Utilities::start_timer( name => "runKrona" );

	#Need to move this to the merge summary function
	unless ($Phylosift::Settings::simple) {

		# skip this if only a simple summary is desired (it's slow)
		debug "Generating krona\n";
		if ($Phylosift::Settings::extended) {
			%ncbi_summary = ();
			%ncbi_summary = %all_summary;
		}
		my $html_report = $Phylosift::Settings::file_dir."/".$self->{"fileName"}.".html";

		#$self->{HTML} = ps_open(">$html_report") unless defined $self->{HTML};
		if ( scalar( keys(%concat_summary) ) > 0 ) {
			if ( !defined $self->{HTML} ) {
				$self->{HTML} = Phylosift::HTMLReport::begin_report( self => $self, file => $html_report );
			}
			Phylosift::HTMLReport::add_krona( self => $self, OUTPUT => $self->{HTML}, summary => \%concat_summary );
		} else {
			my $FH = ps_open(">$html_report");
			print $FH "No results Found\n";
			close($FH);
		}

		#Phylosift::HTMLReport::add_run_info( self => $self, OUTPUT => $self->{HTML}) if exists $self->{HTML};
	}
	Phylosift::Utilities::end_timer( name => "runKrona" );

}

sub print_run_info {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $OUTPUT  = $args{OUTPUT} || miss("OUTPUT");
	my $newline = $args{newline} || "\n";
	print $OUTPUT "Program run as:$newline".join( " ", "phylosift", @{ $self->{"ARGV"} } )."$newline";
	print $OUTPUT "Marker database version:\n"
	  .join( "$newline", Phylosift::Utilities::get_marker_version( path => $Phylosift::Settings::marker_dir ) )
	  ."$newline";
}

#
# non-functional until dependency on Math::Random can be eliminated
#
sub write_confidence_intervals {
	my %args         = @_;
	my $self         = $args{self} || miss("self");
	my $ncbireadsref = $args{ncbi_reads_reference} || miss("ncbi_reads_reference");
	my $totalreads   = $args{total_reads} || miss("total_reads");
	my %ncbireads    = %$ncbireadsref;

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
	for ( my $sI = 0; $sI < $sample_count; $sI++ ) {

		#        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
		my @sample;
		my $kI = 0;
		foreach my $key ( keys(%ncbireads) ) {
			push( @{ $samples{$key} }, $sample[ $kI++ ] );
		}
	}
	my $TAXA_CONF = ps_open( ">".$Phylosift::Settings::file_dir."/taxaconfidence.txt" );
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
		return ( $in, "", "" ) if defined( $deleted->{$in} );
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
	my @t = split( /\s+/, $name );
	my $input_name = join( " ", @t );
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
	if ( !-e $self->{"blastDir"}."/lookup_ID.$chunk.tbl" ) {
		warn "Lookup_ID.$chunk.tbl was not found\n";
		return;
	}

	#read the read mapping file
	my %name_mapping = ();
	my $FH_MAP       = ps_open( $self->{"blastDir"}."/lookup_ID.$chunk.tbl" );
	while (<$FH_MAP>) {
		chomp($_);
		my @line = split( /\t/, $_ );
		$line[1] =~ s/\//./g;
		$name_mapping{ $line[1] } = $line[0];
	}
	close($FH_MAP);

	# when running in isolate mode we want the resulting sequence names to include the filename
	if ($Phylosift::Settings::isolate) {
		foreach my $k ( keys(%name_mapping) ) {
			$name_mapping{$k} = basename( $self->{"readsFile"} ).":".$name_mapping{$k};
		}
	}

	my @array_to_rename = ();
	push( @array_to_rename, glob( $self->{"blastDir"}."/*.ffn.$chunk*" ) );
	push( @array_to_rename, glob( $self->{"blastDir"}."/*.aa.$chunk*" ) );
	push( @array_to_rename, glob( $self->{"alignDir"}."/*.newCandidate.aa.$chunk" ) );
	push( @array_to_rename, glob( $self->{"alignDir"}."/*.$chunk.unmasked" ) );
	push( @array_to_rename, glob( $self->{"alignDir"}."/*.updated.fasta" ) );
	push( @array_to_rename, glob( $self->{"alignDir"}."/*.newCandidate.aa" ) );
	push( @array_to_rename, glob( $self->{"alignDir"}."/*.$chunk.fasta" ) );
	push( @array_to_rename, glob( $self->{"treeDir"}."/*.$chunk.jplace" ) );
	push( @array_to_rename, glob( $Phylosift::Settings::file_dir."/*.jplace" ) );
	push( @array_to_rename, glob( $Phylosift::Settings::file_dir."/sequence_taxa*.$chunk.txt" ) );
	foreach my $file (@array_to_rename) {
		debug "Processing $file\n";
		my $FH  = ps_open($file);
		my $TMP = ps_open(">$file.tmp");
		if ( $file =~ m/\.jplace/ ) {
			my $chunk_size = 100000;	# process the file in chunks
			while (my $bigfish = <$FH>) {
				for(my $i=0; $i < length($bigfish); $i+=$chunk_size){
					my $llama = substr($bigfish, $i, $chunk_size);
					while ( $llama =~ s/"nm":\[\["(\d+\.\d)/"nm":\[\["$name_mapping{$1}/ ) {}
					print $TMP $llama;
				}
			}
		} else {
			while (<$FH>) {
				if(s/^>(\d+\.\d)/>$name_mapping{$1}/){}
				else{
					s/^(\d+\.\d)/$name_mapping{$1}/; #summary files
				}
				print $TMP $_;
			}
		}
		close($FH);
		close($TMP);
		`mv "$file.tmp" "$file"`;
	}
	my $cmd = "rm $self->{\"blastDir\"}"."/lookup_ID.$chunk.tbl";
	`$cmd`;    #remove the lookup table so we don't try to rename sequences after it has already been done.
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
