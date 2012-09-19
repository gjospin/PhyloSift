package Phylosift::pplacer;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(:all);
use Phylosift::Summarize;
use Phylosift::Settings;
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
set_default_values();
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
	if ( defined( $Phylosift::Settings::coverage ) && $Phylosift::Settings::coverage ne ""  ) {
		$covref = Phylosift::Summarize::read_coverage( file => $Phylosift::Settings::coverage );
	}
	unshift(@{$markRef},'concat'); #adds the concatenation to the list of markers
	my $long_rna_switch = -1;
	my $short_rna_jplace;
	for(my $mI=0; $mI < @{$markRef}; $mI++){
		my $marker = @$markRef[$mI];
		# the PMPROK markers are contained in the concat above
		next if($marker =~ /PMPROK/ && $Phylosift::Settings::updated);
		next if($marker =~ /DNGNGWU/ && $Phylosift::Settings::updated);
		# RNA markers have short and long types
		my $protein = Phylosift::Utilities::is_protein_marker(marker=>$marker);
		$long_rna_switch = ($long_rna_switch + 1)%2 unless $protein;
		my $mbname = Phylosift::Utilities::get_marker_basename(marker=>$marker);
		my $short_rna = $protein ? 0 : !$long_rna_switch;
		my $long_rna = $protein ? 0 : $long_rna_switch;
		$mI-- if($short_rna); # need to hit RNA markers twice: short then long
		my $read_alignment_file = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta( marker => $mbname, chunk => $chunk, long => $long_rna, short => $short_rna );
		my $place_file = "";
		if(-e $read_alignment_file && -s $read_alignment_file > 0){
			my $options = $marker eq "concat" ? "--groups $Phylosift::Settings::pplacer_groups" : "";
			$options .= " --mmap-file abracadabra " if ($marker =~ /18s/ || $marker =~ /16s/ || $marker eq "concat");
			$place_file = place_reads(self=>$self, marker=>$marker, dna=>0, chunk => $chunk, reads=>$read_alignment_file, options=>$options, short_rna=>$short_rna);
			unlink("$self->{\"treeDir\"}/abracadabra") if $options =~ /abracadabra/;	# remove the mmap file created by pplacer
		}
		$short_rna_jplace = $place_file if ($short_rna);
		# merge output from short and long RNA
		if($long_rna==1){
			my $dest_jplace = Phylosift::Utilities::get_aligner_output_fasta( marker => $mbname, chunk => $chunk );
			$dest_jplace = $self->{"treeDir"} . "/" . basename( $dest_jplace, ".fasta" );
			$dest_jplace .= ".jplace";
			my $mver;
			$mver = $short_rna_jplace if(-e $short_rna_jplace && !-e $place_file);
			$mver = $place_file if(!-e $short_rna_jplace && -e $place_file);
			`mv $mver $dest_jplace` if defined $mver;
			if(-e $short_rna_jplace && -e $place_file){
				# both short and long RNA seqs, need to merge.
				my $merge_cl = "$Phylosift::Settings::guppy merge -o $dest_jplace $place_file $short_rna_jplace";
				`$merge_cl`;
				unlink($place_file);
				unlink($short_rna_jplace);
			}
		}
	}
	if ( defined($chunk) && $self->{"mode"} eq "all" ) {
		Phylosift::Utilities::end_timer( name => "runPplacer" );
		Phylosift::Utilities::start_timer( name => "runSummarize" );
		Phylosift::Summarize::summarize( self => $self, marker_reference => $markRef , chunk=>$chunk  );
	}
}

=head2 set_default_values

set_default_values for all the parameters in this module

=cut
sub set_default_values{
	my %args = @_;
	my $self = $args{self};
	Phylosift::Settings::set_default(parameter=>\$Phylosift::Settings::pplacer_groups,value=>15);
	Phylosift::Settings::set_default(parameter=>\$Phylosift::Settings::pplacer_verbosity,value=>"0");
	Phylosift::Settings::set_default(parameter=>\$Phylosift::Settings::max_submarker_dist,value=>0.15);
	Phylosift::Settings::set_default(parameter=>\$Phylosift::Settings::min_submarker_prob,value=>0.35);	
}

sub merge_chunk {
	my %args = @_;
	my $chunk = $args{chunk};
	my $place_file = $args{place_file};
	# make sure there's actually work to be done
	return unless defined($chunk) && defined($place_file);
	my $unchunked_place = $place_file;
	$unchunked_place =~ s/\.\d+\.jplace/.jplace/g;
	if(-e $unchunked_place){
		# merge
		my $merge_cl = "$Phylosift::Settings::guppy merge -o $unchunked_place $place_file $unchunked_place";
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
	$edge = $1 if $name=~ /\{(\d+)\}/g;
	return $subgroup_dist{$edge};
}

sub set_submarker_distance {
	my %args  = @_;
	my $name  = $args{name} || miss("name");
	my $dist  = $args{dist};
	my $group = $args{group};
	my $edge  = $name;
	$edge = $1 if $name=~ /\{(\d+)\}/g;
	$subgroup_dist{$edge} = [$dist, $group];
}

my %submarker_map;
sub load_submarkers {
	my %args  = @_;
	return if %submarker_map;
	return unless -e "$Phylosift::Settings::marker_dir/submarkers.txt";
	my $SUBS = ps_open("$Phylosift::Settings::marker_dir/submarkers.txt");
	while(my $line = <$SUBS>){
		chomp $line;
		my @data = split(/\t/,$line);
		$submarker_map{$data[0]}=$data[3];
	}
}

sub get_submarker {
	my %args  = @_;
	my $name  = $args{name} || miss("name");
	$name =~ s/\{\d+\}//g;
	return $submarker_map{$name};
}

sub find_submarker_reads {
	my %args       = @_;
	my $place_file = $args{place_file} || miss("place_file");

	debug "submarker read jplace\n";
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
					my $subinfo = get_submarker_distance(name=>$child->get_name);
					$subinfo->[0] += $child->get_branch_length;
					if( $subinfo->[0] < $min_dist ){
						$min_dist = $subinfo->[0];
						$subgroup = $subinfo->[1];
					}
				}
				set_submarker_distance(name=>$name, dist=>$min_dist, group=>$subgroup);
			}
		}
	);
	
	# submarker criteria:
	# send a read to a submarker if at least X % of its placement probability mass is within distance Y of submarker members
	# 
	my $max_submarker_dist = $Phylosift::Settings::max_submarker_dist;	# rough guesstimate -- could use tuning
	my $min_submarker_prob = $Phylosift::Settings::min_submarker_prob;	# rough guesstimate
	
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
			my $sd = get_submarker_distance( edge=>$edge );
			if($sd->[0] + $pendant < $max_submarker_dist){
				$sub_prob{$sd->[1]} = 0 unless defined($sub_prob{$sd->[1]});
				$sub_prob{$sd->[1]} += $sd->[0] + $pendant;
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
	my $chunk = $args{chunk};

	# determine which reads go to which submarker
	my $subreads = find_submarker_reads(place_file=>$place_file);
	
	# convert to groups
	my %groups;
	foreach my $read (keys %$subreads){
		push( @{$groups{$subreads->{$read}}}, $read );
	}
	
	# filter the codon alignment into subalignments
	my $codon_file = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta( marker => $marker, dna=>1, chunk => $chunk );
	my $alnio = Bio::AlignIO->new(-file => $codon_file );
	my $codon_aln = $alnio->next_aln;
	foreach my $group( keys(%groups)){
		my $sub_file = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta( marker => $marker, sub_marker=>$group, dna=>1, chunk => $chunk );
		my $ALNOUT = ps_open(">$sub_file");
		debug scalar(@{$groups{$group}})." reads in group $group\n";
		debug "reads are ".join(" ",@{$groups{$group}})."\n";
		foreach my $id( @{$groups{$group}}){
			foreach my $seq ( $codon_aln->each_seq_with_id($id) ) {
				print $ALNOUT ">$id\n".$seq->seq."\n";
			}
		}
		close $ALNOUT;
	}
	
	# now place reads on each of these subalignments
	foreach my $group( keys(%groups)){
		my $group_aln = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta( marker => $marker, sub_marker=>$group, dna=>1, chunk => $chunk );
		place_reads(self=>$self, reads=>$group_aln, marker=>$marker, dna=>1, sub_marker=>$group, chunk => $chunk);
	}
}

sub place_reads{
	my %args = @_;
	my $self = $args{self} || miss("self");
	my $dna = $args{dna};
	my $marker = $args{marker} || miss("marker");
	my $reads = $args{reads} || miss("reads");
	my $covref = $args{coverage};
	my $submarker = $args{sub_marker};
	my $chunk = $args{chunk};
	my $options = $args{options} || "";
	my $short_rna = $args{short_rna} || 0;
	my $marker_package = Phylosift::Utilities::get_marker_package( self => $self, marker => $marker, dna => $dna, sub_marker=>$submarker );
	unless(-d $marker_package ){
		# try not updated
		if($Phylosift::Settings::updated){
			$Phylosift::Settings::updated=0;
			$marker_package = Phylosift::Utilities::get_marker_package( self => $self, marker => $marker, dna => $dna, sub_marker=>$submarker );
			$Phylosift::Settings::updated=1;
		}
		unless(-d $marker_package ){
			croak("Marker: $marker\nPackage: $marker_package\nPackage does not exist\nPlacement without a marker package is no longer supported");
		}
	}
	$marker_package .= ".short" if $short_rna;
	my $jplace = $self->{"treeDir"} . "/" . basename($reads, ".fasta").".jplace";
	my $pp = "cd $self->{\"treeDir\"}; $Phylosift::Settings::pplacer $options -o $jplace --verbosity $Phylosift::Settings::pplacer_verbosity -j ".$Phylosift::Settings::threads." -c $marker_package \"$reads\"";
	debug "Running $pp\n";
	system($pp);
	

	return unless -e $jplace;

	# if we're on the concat marker, create a single jplace with all reads for use with multisample metrics 
	if($marker eq "concat"){
		my $sample_jplace = $Phylosift::Settings::file_dir."/".$self->{"fileName"}.".jplace";
		my $sample_jplace_naming = $Phylosift::Settings::file_dir."/".$self->{"fileName"}.".naming.jplace";
		my $markermapfile = Phylosift::Utilities::get_marker_taxon_map(self=>$self, marker=>$marker, dna=>$dna, sub_marker=>$submarker);
		return unless -e $markermapfile;	# can't summarize if there ain't no mappin'!
		my $taxonmap = Phylosift::Summarize::read_taxonmap(file=>$markermapfile);
		name_taxa_in_jplace( self => $self, input => $jplace, output => $sample_jplace_naming, taxonmap=>$taxonmap );
		`$Phylosift::Settings::guppy merge -o $sample_jplace $sample_jplace_naming $sample_jplace` if -f $sample_jplace;
		`cp $sample_jplace_naming $sample_jplace` unless -f $sample_jplace;
		`rm $sample_jplace_naming`;
	}

	if(!$dna && $Phylosift::Settings::updated && Phylosift::Utilities::is_protein_marker(marker=>$marker)){
	    debug "Placing on sub markers $marker\n";
		load_submarkers();
		if(keys(%submarker_map)>0){
			make_submarker_placements(self=>$self, marker=>$marker, chunk=>$chunk, place_file=>$jplace) if $self->{"dna"};
		}
	}
	
	unless($Phylosift::Settings::simple){
		# skip this if a simple summary if desired since it's slow.
		debug "Naming taxa in marker $marker\n";

		# read the tree edge to taxon map for this marker
		my $markermapfile = Phylosift::Utilities::get_marker_taxon_map(self=>$self, marker=>$marker, dna=>$dna, sub_marker=>$submarker);
		return $jplace unless -e $markermapfile;	# can't summarize if there ain't no mappin'!
		my $taxonmap = Phylosift::Summarize::read_taxonmap(file=>$markermapfile);

		# rename nodes
		name_taxa_in_jplace( self => $self, input => $jplace, output => $jplace, taxonmap=>$taxonmap );
	}
	return $jplace unless defined($covref);
	debug "Weighting sequences in $marker\n";
	weight_placements( self => $self, coverage => $covref, place_file => $jplace );
	return $jplace;
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
	`mkdir "$Phylosift::Settings::file_dir"` unless ( -e $Phylosift::Settings::file_dir );
	`mkdir "$self->{"treeDir"}"` unless ( -e $self->{"treeDir"} );
}

=head1 SUBROUTINES/METHODS

=head2 nameTaxa

=cut

my %namemap;
sub read_name_map {
	my %args   = @_;
	return \%namemap if %namemap;
	my $id_file = Phylosift::Utilities::get_gene_id_file();
	my $NAMETABLE = ps_open( $id_file );
	while ( my $line = <$NAMETABLE> ) {
		chomp $line;
		my @dat = split( /\t/, $line );
		$namemap{ $dat[2] } = $dat[1];
	}
	return \%namemap;
}

sub name_taxa_in_jplace {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $input  = $args{input} || miss("input");
	my $output = $args{output} || miss("output");
	my $taxonmap = $args{taxonmap} || miss("taxonmap");

	# read in the taxon name map
	my $namemap = read_name_map();
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
	$tree_string =~ s/:(.+?)(\{\d+?\})/$2:$1/g;
	$tree_string =~ s/\"\,$//g;
	my $tree = Bio::Phylo::IO->parse(
									  '-string' => $tree_string,
									  '-format' => 'newick',
	)->first;

	foreach my $node ( @{ $tree->get_entities } ) {

		my $name = $node->get_name;
		$name =~ s/\{\d+?\}//g;
		if(defined($namemap->{$name})){
			my @data = Phylosift::Summarize::get_taxon_info( taxon => $namemap->{$name} );
			my $ncbi_name = Phylosift::Summarize::tree_name( name => $data[0] );
			$node->set_name($ncbi_name);
			next;
		}
		# this might be an internal node. try to get a taxon group ID and name from the taxon map
		$name = $node->get_name;
		my $edge_id = $1 if $name =~ /\{(\d+?)\}/; 						
		next unless defined($taxonmap->{$edge_id});
		my $node_name="";
		foreach my $tid(@{$taxonmap->{$edge_id}}){
			#debug "TID: $tid\t";
			my @data = Phylosift::Summarize::get_taxon_info( taxon => $tid );
			#debug "$data[0]\n";
			next unless defined $data[0];
			my $ncbi_name = Phylosift::Summarize::tree_name( name => $data[0] ) ;
			$node_name .= "_$ncbi_name"."_" if defined $ncbi_name;
		}
		$node->set_name($node_name);
	}
	$json_data->{tree} = unparse( '-phylo' => $tree, '-format' => 'newick', '-nodelabels' => 1 );
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
