package Phylosift::pplacer;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(:all);
use Phylosift::Summarize;
use Bio::Phylo::IO qw(parse unparse);
use Bio::AlignIO;
use File::Basename;
use JSON;
use Carp;

=head1 NAME

Phylosift::pplacer - place aligned reads onto a phylogenetic tree with pplacer

=head1 VERSION

Version 0.01

=cut
our $VERSION = '0.01';

=head1 SYNOPSIS

Runs Pplacer for all the markers in the list passed as input.
the script will look for the files of interest in
$workingDir/markers/
AND
$workingDir/PS_temp/alignDir/

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 pplacer

=cut

sub pplacer {
	my %args    = @_;
	my $self    = $args{self} || miss("self");
	my $markRef = $args{marker_reference} || miss("marker_reference");
	my $chunk   = $args{chunk};
	
	directoryPrepAndClean( self => $self );
	# if we have a coverage map then weight the placements
	my $covref;
	if ( defined( $self->{"coverage"} ) ) {
		$covref = Phylosift::Summarize::read_coverage( file => $self->{"coverage"} );
	}

	if ( $self->{"updated"} && ! $self->{"extended"} ) {
		my $outputFastaAA = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => "concat", chunk => $chunk );
		my $place_file = place_reads(self=>$self, marker=>"concat", options=>"--groups 10", dna=>0, reads=>$outputFastaAA);
		merge_chunk(chunk => $chunk, place_file=>$place_file);
	}
	foreach my $marker ( @{$markRef} ) {
		# the PMPROK markers are contained in the concat above
		next if($marker =~ /PMPROK/ && $self->{"updated"});
		my $read_alignment_file = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker, chunk => $chunk );
		print "ALIGNMENT FILE : $read_alignment_file\n";
		next unless -e $read_alignment_file;
		
		my $place_file = place_reads(self=>$self, marker=>$marker, dna=>0, reads=>$read_alignment_file);
		# if we're chunked, merge this with the main jplace
		merge_chunk(chunk => $chunk, place_file=>$place_file);
	}
}

sub merge_chunk {
	my %args = @_;
	my $chunk = $args{chunk};
	my $place_file = $args{place_file};
	# make sure there's actually work to be done
	return unless defined($chunk) && defined($place_file);
	my $unchunked_place = $place_file;
	$unchunked_place =~ s/\.\d+\.trim/.trim/g;
	if(-e $unchunked_place){
		# merge
		my $merge_cl = "$Phylosift::Utilities::guppy merge -o $unchunked_place $place_file $unchunked_place";
		debug("Merging jplaces with $merge_cl\n");
		system($merge_cl);
	}else{
		# move
		`mv $place_file $unchunked_place`;
	}
}

my %subgroup_dist;
sub get_submarker_distance {
	my %args = @_;
	my $name = $args{name};
	my $edge = $args{edge};
	return $subgroup_dist{$edge} if defined($edge);
	$edge = $1 if $name=~ /\{\d+\}/g;
	return $subgroup_dist{$edge};
}

sub set_submarker_distance {
	my %args  = @_;
	my $name  = $args{name} || miss("name");
	my $dist  = $args{dist} || miss("dist");
	my $group = $args{group} || miss("group");
	my $edge = $1 if $name=~ /\{\d+\}/g;
	my @dist = ($dist, $group);
	$subgroup_dist{$edge} = @dist;
}

my %submarker_map;
sub load_submarkers {
	my %args  = @_;
	return if %submarker_map;
	my $SUBS = ps_open("submarkers.txt");
	while(my $line = <$SUBS>){
		chomp $line;
		my @data = split(/\t/,$line);
		$submarker_map{$data[0]}=$data[3];
	}
}

sub get_submarker {
	my %args  = @_;
	my $name  = $args{name} || miss("name");
	return $submarker_map{$name};
}

sub find_submarker_reads {
	my %args       = @_;
	my $place_file = $args{place_file} || miss("place_file");

	# read the jplace
	my $JPLACEFILE = ps_open( $place_file );
	my @treedata = <$JPLACEFILE>;	
	close $JPLACEFILE;
	my $json_data = decode_json( join("", @treedata) );
	
	# parse the tree
	my $tree_string = $json_data->{tree};
	$tree_string =~ s/^\s+\"//g;
	$tree_string =~ s/\"\,$//g;
	# move the edge numbers into the node label position
	$tree_string =~ s/:(.+?)\{(\d+?)\}/\{$2\}:$1/g;
	my $tree = Bio::Phylo::IO->parse(
									  '-string' => $tree_string,
									  '-format' => 'newick',
	)->first;
	
	# load the submarkers for this marker
	load_submarkers();
	
	# identify the nearest submarker for each node
	# this is done by depth first search. 
	# at each post-order step the shortest distance to a subgroup member is updated
	my $name_dist_map;
	$tree->visit_depth_first(
		-post => sub {
			my $node = shift;
			my $name = $node->get_name;
			# is this one in a subgroup? 
			# if so set our subgroup distance to 0
			my $subgroup = get_submarker(name=>$name);
			if(defined($subgroup)){
				set_submarker_distance(name=>$name, dist=>0, group=>$subgroup);
			}else{
				my $min_dist = 999999;
				foreach my $child( @{ $node->get_children } ){
					my ($sd, $g) = get_submarker_distance(name=>$child->get_name);
					$sd += $child->get_distance;
					if( $sd < $min_dist ){
						$min_dist = $sd;
						$subgroup = $g;
					}
				}
				set_submarker_distance(name=>$name, dist=>$min_dist, group=>$subgroup);
			}
		}
	);
	
	# submarker criteria:
	# send a read to a submarker if at least X % of its placement probability mass is within distance Y of submarker members
	# 
	my $max_submarker_dist = 0.25;	# rough guesstimate -- could use tuning
	my $min_submarker_prob = 0.5;	# rough guesstimate
	
	# for each placed read, find its probability mass near submarkers and assign to a submarker
	# if appropriate
	my %submarker_reads;
	for(my $i=0; $i< @{$json_data->{placements}}; $i++){
		my $read = $json_data->{placements}->[$i];
		my %sub_prob;
		for( my $j=0; $j < @{$read->{p}}; $j++){
			my $edge = $read->{p}->[$j]->[0];
			my $distal = $read->{p}->[$j]->[3];
			my $pendant = $read->{p}->[$j]->[4];
			my ($sd, $g) = get_submarker_distance( edge=>$edge );
			if($sd + $pendant < $max_submarker_dist){
				$sub_prob{$g} = 0 unless defined($sub_prob{$g});
				$sub_prob{$g} += $sd + $pendant;
			}
		}
		my $submarker;
		foreach my $g(keys(%sub_prob)){
			$submarker = $g if ($sub_prob{$g} > $min_submarker_prob);
		}
		next unless defined($submarker);
		# remove this one from the .jplace
		splice @{$json_data->{placements}}, $i--, 1;

		# add it to the list of reads belonging to submarkers
		for(my $k=0; $k < @{$read->{nm}}; $k++){
			my $qname = $read->{nm}->[$k]->[0];
			$submarker_reads{$qname} = $submarker;
		}
	}	

	# write a new jplace without the submarker reads	
	$JPLACEFILE = ps_open( ">".$place_file );
	print $JPLACEFILE encode_json( $json_data );
	close $JPLACEFILE;
	
	return \%submarker_reads;
}

sub make_submarker_placements{
	my %args       = @_;
	my $self = $args{self} || miss("self");
	my $marker = $args{marker} || miss("marker");
	my $place_file = $args{place_file} || miss("place_file");

	# determine which reads go to which submarker
	my $subreads = find_submarker_reads(place_file=>$place_file);
	
	# convert to groups
	my %groups;
	foreach my $read (keys %$subreads){
		push( @{$groups{$subreads->{$read}}}, $read );
	}
	
	# filter the codon alignment into subalignments
	my $codon_file = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_DNA( marker => $marker );
	my $alnio = Bio::AlignIO->new(-file => $codon_file );
	my $codon_aln = $alnio->next_aln;
	foreach my $group( keys(%groups)){
		my $group_aln = $codon_file. ".sub$group"; # FIXME: this needs to go into a function!!!
		my $ALNOUT = ps_open(">".$group_aln);
		foreach my $id( @{$groups{$group}}){
			foreach my $seq ( $codon_aln->each_seq_with_id($id) ) {
				print $ALNOUT ">$id\n".$seq->seq."\n";
			}
		}
		close $ALNOUT;
	}
	
	# now place reads on each of these subalignments
	foreach my $group( keys(%groups)){
		my $group_aln = $codon_file. ".sub$group"; # FIXME: this needs to go into a function!!!
		place_reads(reads=>$group_aln, marker=>$marker, dna=>1);
	}
}

sub place_reads{
	my %args = @_;
	my $self = $args{self} || miss("self");
	my $dna = $args{dna};
	my $marker = $args{marker} || miss("marker");
	my $reads = $args{reads} || miss("reads");
	my $covref = $args{coverage};
	my $options = $args{options} || "";
	print "Placing for $marker\n";
	my $marker_package = Phylosift::Utilities::get_marker_package( self => $self, marker => $marker );
	print "MARKER PACKAGE : $marker_package\n";
	unless(-d $marker_package ){
		return;
		croak("Marker: $marker\nPackage: $marker_package\nPackage does not exist\nPlacement without a marker package is no longer supported");
	}
	my $pp = "$Phylosift::Utilities::pplacer $options --verbosity 0 -j ".$self->{"threads"}." -c $marker_package \"$reads\"";
	print "Running $pp\n";
	system($pp);
	
	my $jplace = basename($reads, ".fasta").".jplace";

	`mv "$jplace" "$self->{"treeDir"}"` if ( -e $jplace );
	
	return unless -e $self->{"treeDir"} . "/$jplace";
	unless($self->{"simple"}){
		# skip this if a simple summary if desired since it's slow.
		debug "Naming taxa in marker $marker\n";
		name_taxa_in_jplace( self => $self, input => $self->{"treeDir"} . "/$jplace", output => $self->{"treeDir"} . "/$jplace" );
	}
	return $self->{"treeDir"} . "/$jplace" unless defined($covref);
	debug "Weighting sequences in $marker\n";
	weight_placements( self => $self, coverage => $covref, place_file => $self->{"treeDir"} . "/$jplace" );
	return $self->{"treeDir"} . "/$jplace";
}

=head2 weight_placements

=cut

sub weight_placements {
	my %args       = @_;
	my $coverage   = $args{coverage} || miss("coverage");
	my $place_file = $args{place_file} || miss("place_file");

	# read the jplace
	my $JPLACEFILE = ps_open( $place_file );
	my @treedata = <$JPLACEFILE>;	
	close $JPLACEFILE;
	my $json_data = decode_json( join("", @treedata) );

	# re-weight each placement
	for(my $i=0; $i< @{$json_data->{placements}}; $i++){
		my $placement = $json_data->{placements}->[$i];
		for(my $j=0; $j < @{$placement->{nm}}; $j++){
			my $qname = $placement->{nm}->[$j]->[0];
			if(defined($coverage->{$qname})){
				$json_data->{placements}->[$i]->{nm}->[$j]->[1] = $coverage->{$qname};
			}else{
				warn "Unable to find coverage for $qname\n";
			}
		}
	}
	
	# write the weighted jplace
	my $OUTPLACE = ps_open( ">".$place_file );
	print $OUTPLACE encode_json( $json_data );
	close $OUTPLACE;
}

=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean {
	my %args = @_;
	my $self = $args{self} || miss("self");
	
	#create a directory for the Reads file being processed.
	`mkdir "$self->{"fileDir"}"` unless ( -e $self->{"fileDir"} );
	`mkdir "$self->{"treeDir"}"` unless ( -e $self->{"treeDir"} );
}

=head1 SUBROUTINES/METHODS

=head2 nameTaxa

=cut

sub name_taxa_in_jplace {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $input  = $args{input} || miss("input");
	my $output = $args{output} || miss("output");

	# read in the taxon name map
	my %namemap;
	my $taxon_map_file = Phylosift::Utilities::get_marker_taxon_map( self => $self );
	my $NAMETABLE = ps_open( $taxon_map_file );
	while ( my $line = <$NAMETABLE> ) {
		chomp $line;
		my @pair = split( /\t/, $line );
		$namemap{ $pair[0] } = $pair[1];
	}
	Phylosift::Summarize::read_ncbi_taxon_name_map();

	# parse the tree file to get leaf node names
	# replace leaf node names with taxon labels
	my $JPLACEFILE = ps_open( $input );
	my @treedata = <$JPLACEFILE>;	
	close $JPLACEFILE;

	my $json_data = decode_json( join("", @treedata) );
	my $tree_string = $json_data->{tree};
	
	# get rid of the leaf numbers and some other mumbo jumbo
	$tree_string =~ s/^\s+\"//g;
	$tree_string =~ s/\{\d+?\}//g;
	$tree_string =~ s/\"\,$//g;
	my $tree = Bio::Phylo::IO->parse(
									  '-string' => $tree_string,
									  '-format' => 'newick',
	)->first;

	foreach my $node ( @{ $tree->get_entities } ) {

		my $name = $node->get_name;
		next unless defined $namemap{$name};
		my @data = Phylosift::Summarize::get_taxon_info( taxon => $namemap{$name} );
		my $ncbi_name = Phylosift::Summarize::tree_name( name => $data[0] );
		$node->set_name($ncbi_name);
	}
	$json_data->{tree} = unparse( '-phylo' => $tree, '-format' => 'newick' );
	$JPLACEFILE = ps_open( ">$output" );
	print $JPLACEFILE encode_json($json_data);
	close $JPLACEFILE;
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

    perldoc Phylosift::pplacer


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
1;    # End of Phylosift::pplacer.pm
