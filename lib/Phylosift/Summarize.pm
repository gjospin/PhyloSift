package Phylosift::Summarize;
use warnings;
use strict;
use FindBin;
use Phylosift::Utilities qw(:all);
use Carp;
use Bio::Phylo;
use Bio::Phylo::Forest::Tree;
use IO::File;
use XML::Writer;
use JSON;

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
	return ( \%nameidmap, \%idnamemap ) if %nameidmap;
	my $ncbidir = $Phylosift::Utilities::ncbi_dir;
	my $TAXIDS = ps_open( "$ncbidir/names.dmp" );
	while ( my $line = <$TAXIDS> ) {
		chomp $line;
		if ( ( $line =~ /scientific name/ ) || ( $line =~ /synonym/ ) || ( $line =~ /misspelling/ ) ) {
			my @vals = split( /\s+\|\s+/, $line );
			$nameidmap{ homogenize_name_ala_dongying( name => $vals[1] ) } = $vals[0];
			$idnamemap{ $vals[0] } = homogenize_name_ala_dongying( name => $vals[1] ) if ( $line =~ /scientific name/ );
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
	my $ncbidir = $Phylosift::Utilities::ncbi_dir;
	my $TAXSTRUCTURE = ps_open( "$ncbidir/nodes.dmp" );
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
	my $COVERAGE = ps_open( $file );
	while ( my $line = <$COVERAGE> ) {
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
# keep a hash counting up all the read placements
# make this File-scope so anonymous functions below can see it
my %ncbireads;
my %ncbi_summary;
sub summarize {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $markRef = $args{marker_reference} || miss("marker_reference");    # list of the markers we're using
	read_ncbi_taxon_name_map();
	read_ncbi_taxonomy_structure();
	my $markerdir = $Phylosift::Utilities::marker_dir;

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
	my %unclassifiable; # {sequenceID}=mass
	unshift( @{$markRef}, "concat" ) if $self->{"updated"};
	for(my $dna = 0; $dna <2; $dna++){
		
		foreach my $marker ( @{$markRef} ) {
	
			# don't bother with this one if there's no read placements
			my $sub_mark;
			$sub_mark = "*" if $dna;
			my $place_base = $self->{"treeDir"} . "/" . Phylosift::Utilities::get_read_placement_file( marker => $marker, dna => $dna, sub_marker => $sub_mark );
			my @place_files;
			@place_files = glob($place_base) if $dna;	# need to glob on all submarkers if in DNA
			push(@place_files, $place_base) if $dna == 0 && -e $place_base;
			
			foreach my $placeFile(@place_files){
				my $PP_COVFILE = ps_open( ">" . Phylosift::Utilities::get_read_placement_file( marker => $marker ) . ".cov" ) if ( defined $self->{"coverage"} );
	
				my $sub;
				$sub = $1 if $placeFile =~ /\.sub(\d+)\./;
				# first read the taxonomy mapping
				my $markermapfile = Phylosift::Utilities::get_marker_taxon_map(self=>$self, marker=>$marker, dna=>$dna, sub_marker=>$sub);
				next unless -e $markermapfile;	# can't summarize if there ain't no mappin'!
				my $TAXONMAP = ps_open( $markermapfile );
				my %markerncbimap;
				while ( my $line = <$TAXONMAP> ) {
					chomp($line);
					my ( $markerbranch, $ncbiname ) = split( /\t/, $line );
					$markerncbimap{$markerbranch} = [] unless defined( $markerncbimap{$markerbranch} );
					push( @{ $markerncbimap{$markerbranch} }, $ncbiname );
				}
		
				# then read & map the placement
				my $JPLACEFILE = ps_open( $placeFile );
				my @treedata = <$JPLACEFILE>;	
				close $JPLACEFILE;
				my $json_data = decode_json( join("", @treedata) );
		
				# for each placement record
				for(my $i=0; $i< @{$json_data->{placements}}; $i++){
					my $place = $json_data->{placements}->[$i];
					my $qname = $place->{nm}->[0]->[0];
					my $qweight = $place->{nm}->[0]->[1];
					# for each placement edge in the placement record
					for( my $j=0; $j < @{$place->{p}}; $j++){
						my $edge = $place->{p}->[$j]->[0];

						if( !defined($markerncbimap{$edge}) ){
							# mark these reads as unclassifiable
							for(my $k=0; $k < @{$place->{nm}}; $k++){
								my $qname = $place->{nm}->[$k]->[0];
								my $qweight = $place->{nm}->[$k]->[1];
								$unclassifiable{$qname}=0 unless defined($unclassifiable{$qname});
								$unclassifiable{$qname} += $qweight;
							}							
						}

						my $mapcount = scalar( @{ $markerncbimap{$edge} } );
						# for each taxon to which the current phylogeny edge could map
						foreach my $taxon ( @{ $markerncbimap{$edge} } ) {
							my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
							# for each query seq in the current placement record
							for(my $k=0; $k < @{$place->{nm}}; $k++){
								my $qname = $place->{nm}->[$k]->[0];
								my $qweight = $place->{nm}->[$k]->[1];
								$placements{$qname} = () unless defined( $placements{$qname} );
								$placements{$qname}{$taxon_id} = 0 unless defined( $placements{$qname}{$taxon_id} );
								$placements{$qname}{$taxon_id} += $qweight / $mapcount;
								$ncbireads{$taxon} = 0 unless defined $ncbireads{$taxon};
								$ncbireads{$taxon} += $qweight / $mapcount;    # split the p.p. across the possible edge mappings
							}
						}
					}
		
				}
			}
		}

	}

	# make a summary of total reads at each taxonomic level
	# this gets used later in krona output
	foreach my $qname ( keys(%placements) ) {
		my $readsummary = sum_taxon_levels( placements => $placements{$qname} );
		foreach my $taxon_id( keys(%$readsummary) ){
			$ncbi_summary{$taxon_id} = 0 unless defined( $ncbi_summary{$taxon_id} );
			$ncbi_summary{$taxon_id} += $readsummary->{$taxon_id};
		}
	}

	# also write out the taxon assignments for sequences
	my $SEQUENCETAXA = ps_open( ">" . $self->{"fileDir"} . "/sequence_taxa.txt" );
	my $SEQUENCESUMMARY = ps_open( ">" . $self->{"fileDir"} . "/sequence_taxa_summary.txt" );
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
			$taxon_level = "Unknown" unless defined($taxon_level);
			$taxon_name  = "Unknown" unless defined($taxon_name);
			if(exists $self->{"read_names"}{$qname}){
				foreach my $name_ref (@{$self->{"read_names"}{$qname}}){
					#my @name_array = @{$name_ref};
					print $SEQUENCETAXA "$name_ref\t$taxon_id\t$taxon_level\t$taxon_name\t" . $placements{$qname}{$taxon_id} . "\n";
				}
			}
		}
		my $readsummary = sum_taxon_levels( placements => $placements{$qname} );
		foreach my $taxon_id ( sort { $readsummary->{$b} <=> $readsummary->{$a} } keys %{$readsummary} ) {
			my ( $taxon_name, $taxon_level, $tid ) = get_taxon_info( taxon => $taxon_id );
			$taxon_level = "Unknown" unless defined($taxon_level);
			$taxon_name  = "Unknown" unless defined($taxon_name);
			if(exists $self->{"read_names"}{$qname}){
				foreach my $name_ref (@{$self->{"read_names"}{$qname}}){
					#my @name_array = @{$name_ref};
					print $SEQUENCESUMMARY "$name_ref\t$taxon_id\t$taxon_level\t$taxon_name\t" . $readsummary->{$taxon_id} . "\n";
				}
			}
			
		}
	}
	close($SEQUENCESUMMARY);
	close($SEQUENCETAXA);

	# sort descending
	my $TAXAOUT = ps_open( ">" . $self->{"fileDir"} . "/taxasummary.txt" );
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
	debug "Total reads placed is ".scalar(keys(%placements))."\n";
	debug "Total classifiable probability mass is $totalreads\n";

	# write the taxa with 90% highest posterior density, assuming each read is an independent observation
	my $taxasum = 0;
	my $TAXAHPDOUT = ps_open( ">" . $self->{"fileDir"} . "/taxa_90pct_HPD.txt" );
	foreach my $taxon ( sort { $ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads ) {
		$taxasum += $ncbireads{$taxon};
		my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $taxon );
		print $TAXAHPDOUT join( "\t", $taxon_id, $taxon_level, $taxon_name, $ncbireads{$taxon} ), "\n";
		last if $taxasum >= $totalreads * 0.9;
	}
	close($TAXAHPDOUT);
	
	unless($self->{"simple"}){
		# skip this if only a simple summary is desired (it's slow)
		debug "Generating krona\n";
		krona_report(self=>$self);
	}
}

my $xml;
my $KRONA_THRESHOLD = 0.01;
sub krona_report {
	my %args = @_;
	my $self = $args{self} || miss("self");
	
	my $OUTPUT = IO::File->new(">".$self->{"fileDir"}."/krona.html");
	print $OUTPUT <<EOF
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
 <head>
  <meta charset="utf-8"/>
  <base href="http://krona.sourceforge.net/" target="_blank"/>
  <link rel="shortcut icon" href="img/favicon.ico"/>
  <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"http://krona.sourceforge.net\"."}</script>
  <script src="src/krona-2.0.js"></script>
 </head>
 <body>
  <img id="hiddenImage" src="img/hidden.png" style="display:none"/>
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style="display:none">	

EOF
;
	debug "init xml\n";
	$xml = new XML::Writer(OUTPUT=>$OUTPUT);
	$xml->startTag("krona", "collapse"=>"false", "key"=>"true");
	$xml->startTag("attributes", "magnitude"=>"abundance");
	$xml->startTag("attribute", "display"=>"reads");
	$xml->characters("abundance");
	$xml->endTag("attribute");
	$xml->endTag("attributes");
	$xml->startTag("datasets");
	$xml->startTag("dataset");
	$xml->characters($self->{"readsFile"});
	$xml->endTag("dataset");
	$xml->endTag("datasets");

	debug "parse ncbi\n";
	# FIXME: work with other taxonomy trees
	my $taxonomy = Bio::Phylo::IO->parse(
									  '-file'   => "$Phylosift::Utilities::marker_dir/ncbi_tree.updated.tre",
									  '-format' => 'newick',
	)->first;
	
	debug "visitor\n";
	my $root = $taxonomy->get_root;
	debug "Root node id ".$root->get_name."\n";
	debug "Root node read count ".$ncbi_summary{$root->get_name}."\n";
	# write out abundance for nodes that have > $KRONA_THRESHOLD probability mass
	$root->visit_depth_first(
		-pre => sub {
			my $node = shift;
			my $name = $node->get_name;
			return unless(defined($ncbi_summary{$name}) && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD);
			my ( $taxon_name, $taxon_level, $taxon_id ) = get_taxon_info( taxon => $name );
			$xml->startTag("node", "name"=>$taxon_name, "href"=>"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$name");
			$xml->startTag("abundance");
			$xml->startTag("val");
			$xml->characters($ncbi_summary{$name});
			$xml->endTag("val");
			$xml->endTag("abundance");
		},
		-pre_sister => sub {
			my $node = shift;
			my $name = $node->get_name;
			return unless(defined($ncbi_summary{$name}) && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD);
			$xml->endTag("node");
		},
		-no_sister => sub {
			my $node = shift;
			my $name = $node->get_name;
			return unless(defined($ncbi_summary{$name}) && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD);
			$xml->endTag("node");
		}
	);
	debug "done visiting!\n";
	
	$xml->endTag("krona");
	$xml->end();
	print $OUTPUT "\n</div>\n";
	print_run_info(self=>$self, OUTPUT=>$OUTPUT, newline=>"<br/>\n");
	print $OUTPUT "</body></html>\n";
	$OUTPUT->close();
}

sub print_run_info{
	my %args = @_;
	my $self = $args{self} || miss("self");
	my $OUTPUT = $args{OUTPUT} || miss("OUTPUT");
	my $newline = $args{NEWLINE} || "\n";
	print $OUTPUT "Program run as:$newline".join(" ", "phylosift", @{ $self->{"ARGV"} })."$newline";
	print $OUTPUT "Marker database version:\n".join("$newline", Phylosift::Utilities::get_marker_version(path=>$Phylosift::Utilities::marker_dir))."$newline";
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
	for ( my $sI = 0 ; $sI < $sample_count ; $sI++ ) {

		#        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
		my @sample;
		my $kI = 0;
		foreach my $key ( keys(%ncbireads) ) {
			push( @{ $samples{$key} }, $sample[ $kI++ ] );
		}
	}
	my $TAXA_CONF = ps_open( ">" . $self->{"fileDir"} . "/taxaconfidence.txt" );
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
			last if defined($parent{$cur_tid}[0]) && $parent{$cur_tid}[0] == $cur_tid;
			$cur_tid = $parent{$cur_tid}[0];
		}
	}
	return \%summarized;
}

sub get_taxon_info {
	my %args = @_;
	my $in   = $args{taxon} || miss("taxon");
	read_ncbi_taxonomy_structure();
	if ( $in =~ /^\d+$/ ) {

		#it's an ncbi taxon id.  look up its name and level.
		my $merged = Phylosift::Summarize::read_merged_nodes();
		$in = $merged->{$in} if defined($merged->{$in});
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
	my %args   = @_;
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
	my $MERGED = ps_open( "$Phylosift::Utilities::ncbi_dir/merged.dmp" );
	while ( my $line = <$MERGED> ) {
		chomp $line;
		my @vals = split( /\s+\|\s*/, $line );
		$merged{ $vals[0] } = $vals[1];
	}
	debug "Done reading merged\n";
	return \%merged;
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
