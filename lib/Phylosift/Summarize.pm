package Phylosift::Summarize;
use warnings;
use strict;
use FindBin;
use Phylosift::Utilities qw(debug);
use Carp;
use Bio::Phylo;
use Bio::Phylo::Forest::Tree;

if ( $^O =~ /arwin/ ) {
	use lib "$FindBin::Bin/osx/darwin-thread-multi-2level/";
}

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

=head2 read_ncbi_taxon_name_map

# read the NCBI taxon names
# stash them in hashes called nameidmap and idnamemap to go back & forth from tax ids to names

=cut

sub read_ncbi_taxon_name_map {
	return if %nameidmap;
	my $ncbidir = $Phylosift::Utilities::ncbi_dir;
	open( TAXIDS, "$ncbidir/names.dmp" );
	while ( my $line = <TAXIDS> ) {
		chomp $line;
		if ( ( $line =~ /scientific name/ ) || ( $line =~ /synonym/ ) || ( $line =~ /misspelling/ ) ) {
			my @vals = split( /\s+\|\s+/, $line );
			$nameidmap{ homogenize_name_ala_dongying( name => $vals[1] ) } = $vals[0];
			$idnamemap{ $vals[0] } = homogenize_name_ala_dongying( name => $vals[1] ) if ( $line =~ /scientific name/ );
		}
	}
	return ( %nameidmap, %idnamemap );
}

# now read the NCBI taxonomy structure
# puts the results in a hash called "parent"
my %parent;

sub read_ncbi_taxonomy_structure {
	my $ncbidir = $Phylosift::Utilities::ncbi_dir;
	open( TAXSTRUCTURE, "$ncbidir/nodes.dmp" ) || croak "Unable to read $ncbidir/nodes.dmp";
	while ( my $line = <TAXSTRUCTURE> ) {
		chomp $line;
		my @vals = split( /\s+\|\s+/, $line );
		$parent{ $vals[0] } = [ $vals[1], $vals[2] ];
	}
	return %parent;
}

sub make_ncbi_tree_from_update {
	my %args        = @_;
	my $self        = $args{self} // miss("self");
	my $results_dir = $args{results_directory} // miss("results_directory");
	my $markerdir   = $args{marker_directory} // miss("marker_directory");
	read_ncbi_taxon_name_map();
	read_ncbi_taxonomy_structure();
	open( AAIDS,          "$markerdir/gene_ids.aa.txt" );
	open( MARKERTAXONMAP, ">$markerdir/marker_taxon_map.updated.txt" );
	my @taxonids;

	while ( my $line = <AAIDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		push( @taxonids, $taxon ) if $taxon =~ /^\d+$/;
		print MARKERTAXONMAP "$uniqueid\t$taxon\n";
	}
	close MARKERTAXONMAP;
	my %tidnodes;
	my $phylotree = Bio::Phylo::Forest::Tree->new();
	foreach my $tid (@taxonids) {
		next if ( $tid eq "" );
		my $child;
		while ( $tid != 1 ) {

			# check if we've already seen this one
			last if ( defined( $tidnodes{$tid} ) );

			# create a new node & add to tree
			my $parentid = $parent{$tid}->[0];
			my $newnode;
			$newnode = Bio::Phylo::Forest::Node->new( -parent => $tidnodes{$parentid}, -name => $tid ) if defined( $tidnodes{$parentid} );
			$newnode = Bio::Phylo::Forest::Node->new( -name => $tid ) if !defined( $tidnodes{$parentid} );
			$tidnodes{$tid} = $newnode;
			$newnode->set_child($child) if ( defined($child) );
			$phylotree->insert($newnode);

			# continue traversal toward root
			$tid   = $parentid;
			$child = $newnode;
		}
	}
	open( TREEOUT, ">ncbi_tree.updated.tre" );
	print TREEOUT $phylotree->to_newick( "-nodelabels" => 1 );
	close TREEOUT;
}

=head2 make_ncbi_tree
Reads all the marker gene trees, finds their corresponding taxa in the NCBI taxonomy, and
constructs a newick format tree representation of the NCBI taxonomy containing only the
organisms present in the marker gene trees.  
=cut

sub make_ncbi_tree {
	my %args = @_;
	my $self = $args{self} // miss("self");
	read_ncbi_taxon_name_map();
	read_ncbi_taxonomy_structure();

	# now read the list of organisms we have in our DB
	# construct a phylo tree with the NCBI topology containing
	# just the organisms in our database
	my $markerdir = $Phylosift::Utilities::marker_dir;
	my %namemap   = Phylosift::Utilities::read_name_table( marker_directory => $markerdir );
	my $phylotree = Bio::Phylo::Forest::Tree->new();
	open( MARKERTAXONMAP, ">$markerdir/marker_taxon_map.txt" );
	my %tidnodes;
	foreach my $key ( keys(%namemap) ) {
		$namemap{$key} = homogenize_name_ala_dongying( name => $namemap{$key} );
		my ( $tid, $name ) = donying_find_name_in_taxa_db( name => $namemap{$key} );
		if ( $tid eq "ERROR" ) {
			print STDERR "Error! Could not find $namemap{$key} in name map\n" if length($key) > 12;
			next;
		}

		# add it to the mapping file
		my $treename = tree_name( name => $idnamemap{$tid} );
		print MARKERTAXONMAP "$key\t$treename\n";

		#got the taxon id, now walk to root adding tree nodes as necessary
		next unless ( defined($tid) );
		my $child;
		while ( $tid != 1 ) {

			# check if we've already seen this one
			last if ( defined( $tidnodes{$tid} ) );

			# create a new node & add to tree
			my $nodename = tree_name( name => $idnamemap{$tid} );
			my $parentid = $parent{$tid}->[0];
			my $newnode;
			$newnode = Bio::Phylo::Forest::Node->new( -parent => $tidnodes{$parentid}, -name => $nodename ) if defined( $tidnodes{$parentid} );
			$newnode = Bio::Phylo::Forest::Node->new( -name => $nodename ) if !defined( $tidnodes{$parentid} );
			$tidnodes{$tid} = $newnode;
			$newnode->set_child($child) if ( defined($child) );
			$phylotree->insert($newnode);

			# continue traversal toward root
			$tid   = $parentid;
			$child = $newnode;
		}
	}
	close MARKERTAXONMAP;
	open( TREEOUT, ">ncbi_tree.tre" );
	print TREEOUT $phylotree->to_newick( "-nodelabels" => 1 );
	close TREEOUT;
}

=head2 read_coverage
Reads a coverage file
Input: file - a file name
=cut

sub read_coverage {
	my %args = @_;
	my $file = $args{file} // miss("file");
	my %coverage;
	open( COVERAGE, $file ) || croak("Unable to read coverage file $file");
	while ( my $line = <COVERAGE> ) {
		chomp $line;
		my @data = split( /\t/, $line );
		$data[0] =~ s/[\(\)\[\]\+\=\<\>\?]//g;
		$data[0] =~ s/[\-\s]/_/g;
		$coverage{ $data[0] } = $data[1];
	}
	return \%coverage;
}

=head2 summarize
Reads the .place files containing Pplacer read placements and maps them onto the
NCBI taxonomy
=cut

sub summarize {
	my %args    = @_;
	my $self    = $args{self} // miss("self");
	my $markRef = $args{marker_reference} // miss("marker_reference");    # list of the markers we're using
	read_ncbi_taxon_name_map();
	read_ncbi_taxonomy_structure();
	my $markerdir = $Phylosift::Utilities::marker_dir;
	my %namemap = Phylosift::Utilities::read_name_table( marker_directory => $markerdir );
	foreach my $key ( keys(%namemap) ) {
		$namemap{$key} = homogenize_name_ala_dongying( name => $namemap{$key} );
	}

	# keep a hash counting up all the read placements
	my %ncbireads;

	# try to read a contig coverage file if it exists
	my %coverage;
	if ( defined $self->{"coverage"} ) {
		my $covref = read_coverage( file => $self->{"coverage"} );
		%coverage = %$covref;
	}

	# read all of the .place files for markers
	# map them onto the ncbi taxonomy
	# this is a hash structured as {sequenceID}{taxonID}=probabilitySum
	my %placements;
	unshift( @{$markRef}, "concat" ) if $self->{"updated"};
	foreach my $marker ( @{$markRef} ) {
		my $markermapfile = "$markerdir/$marker.ncbimap";
		$markermapfile = "$markerdir/$marker.updated/$marker.taxonmap" if $self->{"updated"};
		next unless -e $markermapfile;

		# don't bother with this one if there's no read placements
		my $placeFile = $self->{"treeDir"} . "/" . Phylosift::Utilities::get_read_placement_file( marker => $marker );
		next unless ( -e $placeFile );
		my $pp_covfile;
		open( $pp_covfile, ">" . Phylosift::Utilities::get_read_placement_file( marker => $marker ) . ".cov" ) if ( defined $self->{"coverage"} );

		# first read the taxonomy mapping
		open( TAXONMAP, $markermapfile ) || croak("Unable to read file $markermapfile\n");
		my %markerncbimap;
		while ( my $line = <TAXONMAP> ) {
			chomp($line);
			my ( $markerbranch, $ncbiname ) = split( /\t/, $line );
			$markerncbimap{$markerbranch} = [] unless defined( $markerncbimap{$markerbranch} );
			push( @{ $markerncbimap{$markerbranch} }, $ncbiname );
		}

		# then read & map the placement
		open( PLACEFILE, $placeFile ) || croak("Unable to read file $placeFile\n");
		my $placeline = 0;
		my %curplaces;
		while ( my $line = <PLACEFILE> ) {
			print $pp_covfile $line if ( defined $self->{"coverage"} );
			$placeline = 1 if ( $line =~ /"placements"/ );
			next if ( $line =~ /^\>/ );
			next if ( $line =~ /^\s*\#/ );
			if ( $placeline == 1 && $line =~ /\[(\d+),\s.\d+\.?\d+,\s(\d+\.?\d*),/ ) {
				$curplaces{$1} = $2;
			}

			# have we reached the end of a placement entry?
			if ( $placeline == 1 && $line =~ /\"nm\"\:\s+\[\[\"(.+?)\",\s+\d+\]\]/ ) {
				my $qname = $1;

				# tally up the probabilities on different NCBI groups
				foreach my $edge ( keys(%curplaces) ) {
					my $weightRatio = $curplaces{$edge};
					$weightRatio *= $coverage{$qname} if defined( $coverage{$qname} );
					print STDERR "Marker $marker missing mapping from phylogeny edge $edge to taxonomy\n" unless defined( $markerncbimap{$edge} );
					next unless defined( $markerncbimap{$edge} );
					my $mapcount = scalar( @{ $markerncbimap{$edge} } );
					print STDERR "Found 0 taxa for edge $edge\n" if ( scalar( @{ $markerncbimap{$edge} } ) == 0 );
					foreach my $taxon ( @{ $markerncbimap{$edge} } ) {
						my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
						$placements{$qname} = () unless defined( $placements{$qname} );
						$placements{$qname}{$taxon_id} = 0 unless defined( $placements{$qname}{$taxon_id} );
						$placements{$qname}{$taxon_id} += $curplaces{$edge} / $mapcount;
						$ncbireads{$taxon} = 0 unless defined $ncbireads{$taxon};
						$ncbireads{$taxon} += $weightRatio / $mapcount;    # split the p.p. across the possible edge mappings
					}

					# if we have read coverage information, add it to an updated placement file
					my $mass = 1;
					$mass = $coverage{$qname} if defined( $coverage{$qname} );
					my $massline = ", \"m\": \"$mass\"\n";
					print $pp_covfile $line if ( defined $self->{"coverage"} );
				}
				%curplaces = ();
			}
		}
		close $pp_covfile if ( defined $self->{"coverage"} );
	}

	# also write out the taxon assignments for sequences
	open( SEQUENCETAXA,    ">" . $self->{"fileDir"} . "/sequence_taxa.txt" );
	open( SEQUENCESUMMARY, ">" . $self->{"fileDir"} . "/sequence_taxa_summary.txt" );
	foreach my $qname ( keys(%placements) ) {

		# sum up all placements for this sequence, use to normalize
		my $placecount = 0;
		foreach my $taxon_id ( keys( %{ $placements{$qname} } ) ) {
			$placecount += $placements{$qname}{$taxon_id};
		}

		# normalize to probability distribution
		foreach my $taxon_id ( sort { $placements{$qname}{$b} <=> $placements{$qname}{$a} } keys %{ $placements{$qname} } ) {
			$placements{$qname}{$taxon_id} /= $placecount;
			my ( $taxon_name, $taxon_level, $tid ) = get_taxon_info( taxon => $taxon_id );
			print SEQUENCETAXA "$qname\t$taxon_id\t$taxon_level\t$taxon_name\t" . $placements{$qname}{$taxon_id} . "\n";
		}
		my $readsummary = sum_taxon_levels( placements => $placements{$qname} );
		foreach my $taxon_id ( sort { $readsummary->{$b} <=> $readsummary->{$a} } keys %{$readsummary} ) {
			my ( $taxon_name, $taxon_level, $tid ) = get_taxon_info( taxon => $taxon_id );
			print SEQUENCESUMMARY "$qname\t$taxon_id\t$taxon_level\t$taxon_name\t" . $readsummary->{$taxon_id} . "\n";
		}
	}
	close(SEQUENCESUMMARY);
	close(SEQUENCETAXA);

	# sort descending
	open( taxaOUT, ">" . $self->{"fileDir"} . "/taxasummary.txt" );
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		$taxon_level = "Unknown" unless defined($taxon_level);
		$taxon_name  = "Unknown" unless defined($taxon_name);
		print taxaOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
	}
	close(taxaOUT);

	# sample from multinomial to get confidence limits
	# get total read count
	my $totalreads = 0;
	foreach my $val ( values(%ncbireads) ) {
		$totalreads += $val;
	}
	debug "Total reads are $totalreads\n";

	# write the taxa with 90% highest posterior density, assuming each read is an independent observation
	my $taxasum = 0;
	open( TAXAHPDOUT, ">" . $self->{"fileDir"} . "/taxa_90pct_HPD.txt" );
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		$taxasum += $ncbireads{$taxon};
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		print TAXAHPDOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
		last if $taxasum >= $totalreads * 0.9;
	}
	close(TAXAHPDOUT);
}

#
# non-functional until dependency on Math::Random can be eliminated
#
sub write_confidence_intervals {
	my %args         = @_;
	my $self         = $args{self} // miss("self");
	my $ncbireadsref = $args{ncbi_reads_reference} // miss("ncbi_reads_reference");
	my $totalreads   = $args{total_reads} // miss("total_reads");
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
	for ( my $sI = 0 ; $sI < $sample_count ; $sI++ ) {

		#        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
		my @sample;
		my $kI = 0;
		foreach my $key ( keys(%ncbireads) ) {
			push( @{ $samples{$key} }, $sample[ $kI++ ] );
		}
	}
	open( taxaCONF, ">" . $self->{"fileDir"} . "/taxaconfidence.txt" );
	foreach my $key ( keys(%samples) ) {
		my @svals = @{ $samples{$key} };
		my @sorted = sort { $a <=> $b } @svals;
		my ( $taxon_name, $taxon_level, $taxon_id ) = getTaxonInfo($key);
		print taxaCONF join( "\t",
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
	my $placements = $args{placements} // miss("placements");
	my %summarized = ();
	foreach my $taxon_id ( keys %$placements ) {
		my $cur_tid = $taxon_id;
		while ( defined($cur_tid) && $cur_tid != 1 ) {
			$summarized{$cur_tid} = 0 unless defined( $summarized{$cur_tid} );
			$summarized{$cur_tid} += $placements->{$taxon_id};
			$cur_tid = $parent{$cur_tid}[0];
		}
	}
	return \%summarized;
}

sub get_taxon_info {
	my %args = @_;
	my $in   = $args{taxon} // miss("taxon");
	if ( $in =~ /^\d+$/ ) {

		#it's an ncbi taxon id.  look up its name and level.
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
	} else {
	}
	return ( $in, "", "" );
}

sub tree_name {
	my %args   = @_;
	my $inName = $args{name} // miss("name");
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
	my %args   = @_;
	my $inName = $args{name} // miss("name");
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
	my $name = $args{name} // miss("name");
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
