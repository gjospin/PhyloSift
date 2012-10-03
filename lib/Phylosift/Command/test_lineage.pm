package Phylosift::Command::test_lineage;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::Simulations;
use Carp;
use Phylosift::Utilities qw(debug);

sub description {
	return "phylosift test_lineage - conduct a statistical test (a Bayes factor) for the presence of a particular lineage in a sample";
}

sub abstract {
	return "conduct a statistical test (a Bayes factor) for the presence of a particular lineage in a sample";
}

sub usage_desc { "test_lineage %o" }

sub options {
	return (
		[ "sample=s",  "Path to a directory containing a phylosift analysis of a sample", { required => 1 }],
		[ "taxon=s@",  "Taxon ID of lineage to test. If more than one taxon is specified, the test applies to all lineages under their most recent common ancestor", { required => 1 }],
	);
}

sub validate {
	my ($self, $opt, $args) = @_;
	
	$self->usage_error("benchmark requires a phylosift output directory") unless @$args > 0;
}

sub load_opt {
	my %args = @_;
	my $opt = $args{opt};
	$Phylosift::Settings::configuration = $opt->{config};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::my_debug = $opt->{debug};

	$Phylosift::Utilities::debuglevel = $Phylosift::Settings::my_debug || 0;	
}

sub execute {
	my ($self, $opt, $args) = @_;
	load_opt(opt=>$opt);
	Phylosift::Command::sanity_check();
	my $ps = new Phylosift::Phylosift();
	Phylosift::Utilities::program_checks();
	Phylosift::Utilities::data_checks( self => $ps );
	
	# convert taxon IDs to a hash
	my %taxon_ids;
	foreach my $tid(@$opt->{taxon}){
		$taxon_ids{$tid}=1;
	}
	
	# get taxon names
	my $namemap = read_name_map();
	Phylosift::Summarize::read_ncbi_taxon_name_map();

	# read in the JSON file
	my $sample_jplace = $opt->{sample}."/treeDir/concat.1.jplace";
	my $JPLACEFILE = ps_open( $sample_jplace );
	my @treedata = <$JPLACEFILE>;	
	close $JPLACEFILE;

	my $json_data = decode_json( join("", @treedata) );
	my $tree_string = $json_data->{tree};

	# find the nodes in question
	my $tree = Bio::Phylo::IO->parse( '-string' => $tree_string, '-format' => 'newick')->first;
	my @nodes;
	foreach my $node ( @{ $tree->get_entities } ) {
		my $name = $node->get_name;
		$name =~ s/\{\d+?\}//g;

		if(defined($namemap->{$name})){
			push (@nodes, $node) if defined( $taxon_ids{ $namemap->{$name} } );
			next;
		}
	}
	
	# get the mrca of nodes in question
	my $mrca = $tree->get_mrca(\@nodes);
	
	# accumulate list of all subtree nodes of interest
	my %subtree_edges;
	$mrca->visit_depth_first(
		-post => sub {
			my $node = shift;
			my $name = $node->get_name;
			$subtree_edges{$1} = 1 if $name =~ /\{(\d+?)\}/;			
		}
	);

	# now walk the list of placements and add up probability mass
	my $bf_numer = 1;
	for ( my $i = 0 ; $i < @{ $json_data->{placements} } ; $i++ ) {
		my $place = $json_data->{placements}->[$i];

		# for each placement edge in the placement record
		my $mass = 0;
		for ( my $j = 0 ; $j < @{ $place->{p} } ; $j++ ) {
			my $edge      = $place->{p}->[$j]->[0];
			next unless defined($subtree_edges{$edge});
			my $mass += $place->{p}->[$j]->[2];
		}
		$bf_numer *= (1.0-$mass);
	}
	print "Hypothesis: taxa in this group have zero abundance\n";
	print "Bayes factor: ". ($bf_numer / (1-$bf_numer)). "\n";
	print "1-3\tBarely worth mentioning\n";
	print "3-10\tSubstantial\n";
	print "10-30\tStrong\n";
	print "30-100\tVery Strong\n";
	print ">100\tDecisive\n";
}

1;