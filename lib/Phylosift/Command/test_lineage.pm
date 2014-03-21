package Phylosift::Command::test_lineage;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use JSON;
use Carp;
use Phylosift::Utilities qw(debug ps_open);

our $VERSION = "v1.0.1";
sub description {
	return "phylosift test_lineage - conduct a statistical test (a Bayes factor) for the presence of a particular lineage in a sample";
}

sub abstract {
	return "conduct a statistical test (a Bayes factor) for the presence of a particular lineage in a sample";
}

sub usage_desc { "test_lineage %o" }

sub options {
	return (
			 [ "codon", "Conduct the lineage test on a codon-based tree. Can give higher resolution in well sampled clades." ],
			 [ "sample=s", "Path to a directory containing a phylosift analysis of a sample", { required => 1 } ],
			 [
				"taxon=s@",
				"Taxon ID of lineage to test. If more than one taxon is specified, the test applies to all lineages under their most recent common ancestor",
				{ required => 1 }
			 ],
			 [ "marker=s", "Apply test to the specified marker gene family", { default => "concat" } ],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
}

sub load_opt {
	my %args = @_;
	my $opt  = $args{opt};
	$Phylosift::Settings::configuration        = $opt->{config};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::my_debug             = $opt->{debug};

	$Phylosift::Utilities::debuglevel = $Phylosift::Settings::my_debug || 0;
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	load_opt( opt => $opt );
	Phylosift::Command::sanity_check();
	my $ps = new Phylosift::Phylosift();
	Phylosift::Utilities::program_checks();
	Phylosift::Utilities::data_checks( self => $ps );

	# convert taxon IDs to a hash
	my %taxon_ids;
	foreach my $tid ( @{ $opt->{taxon} } ) {
		$taxon_ids{$tid} = 1;
	}

	# get taxon names
	Phylosift::Summarize::read_ncbi_taxon_name_map();

	# read in the JSON file
	my $jplace =
	  Phylosift::Utilities::get_read_placement_file( marker => $opt->{marker}, chunk => 1, updated => 1, dna => $opt->{codon}, sub_marker => $opt->{codon} );
	my $sample_jplace = $opt->{sample}."/treeDir/$jplace";
	my $JPLACEFILE    = ps_open($sample_jplace);
	my @treedata      = <$JPLACEFILE>;
	close $JPLACEFILE;

	my $json_data = decode_json( join( "", @treedata ) );
	my $tree_string = $json_data->{tree};

	# find the nodes in question
	my $tree = Bio::Phylo::IO->parse( '-string' => $tree_string, '-format' => 'newick' )->first;
	my @nodes;
	my %subtree_edges;
	foreach my $node ( @{ $tree->get_entities } ) {
		my $name = $node->get_name;
		my $tid = $1 if $name =~ /\[(\d+?)\]\{/;
		$subtree_edges{$1} = 1 if defined( $taxon_ids{$tid} ) && $name =~ /\{(\d+?)\}/;
		push( @nodes, $node ) if defined( $taxon_ids{$tid} );
	}
	debug "Found ".scalar(@nodes)." nodes\n";
	if ( @nodes == 0 ) {
		die( "Error, requested taxa not found in gene family ".$opt->{marker}."\nUnable to apply Bayes factor test\n\n" );
	}

	# get the mrca of nodes in question
	my @mrca_nodes = @nodes;
	my $mrca       = $tree->get_mrca( \@mrca_nodes );
	debug "nodes[0] name is ".$nodes[0]->get_name."\n";
	debug "mrca name is ".$mrca->get_name."\n";

	# if nontrivial mrca then
	# accumulate a list of all subtree nodes of interest
	if ( scalar(@nodes) > 1 ) {
		$mrca->visit_depth_first(
			-post => sub {
				my $node = shift;
				my $name = $node->get_name;
				$subtree_edges{$1} = 1 if $name =~ /\{(\d+?)\}/;
			}
		);
	}
	debug "After including mrca subtree edges we have ".scalar( keys(%subtree_edges) )." target(s)\n";
	debug "Looking on edges ".join( "\t", keys(%subtree_edges) )."\n";

	# now walk the list of placements and add up probability mass
	my $bf_numer = 1.0;
	for ( my $i = 0; $i < @{ $json_data->{placements} }; $i++ ) {
		my $place = $json_data->{placements}->[$i];

		# for each placement edge in the placement record
		my $mass = 0;
		for ( my $j = 0; $j < @{ $place->{p} }; $j++ ) {
			my $edge = $place->{p}->[$j]->[1];
			die
			  "Error, test_lineage requires posterior probabilities of branch placement, please use phylosift to reanalyze the data with the --bayes option.\n\n"
			  if scalar( @{ $place->{p}->[$j] } ) < 6;
			next unless defined( $subtree_edges{$edge} );
			$mass += $place->{p}->[$j]->[5];
		}
		$bf_numer *= ( 1.0 - $mass );
	}

	print "\nNull hypothesis: taxa in this group have zero abundance in the sample\n";
	my $bf = "Infinite -- at least one sequence has 100% probability of deriving from the test lineage(s)";
	$bf = $bf_numer == 1 ? "Zero -- target is beyond limit of detection" : $bf;
	$bf = ( $bf_numer != 0 && $bf_numer != 1 ) ? ( ($bf_numer) / ( 1 - $bf_numer ) ) : $bf;
	print "\nBayes factor: $bf\n\n";
	print "Strength of null hypothesis rejection:\n";
	print "< 1\t Null hypothesis supported\n";
	print "1-3\tBarely worth mentioning\n";
	print "3-10\tSubstantial\n";
	print "10-30\tStrong\n";
	print "30-100\tVery Strong\n";
	print ">100\tDecisive\n\n";
}

1;
