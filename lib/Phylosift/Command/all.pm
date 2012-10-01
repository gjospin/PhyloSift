package Phylosift::Command::all;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Carp;
use Phylosift::Utilities qw(debug);

sub description {
	return "phylosift all - run all steps for phylogenetic analysis of genomic or metagenomic sequence data";
}

sub abstract {
	return "run all steps for phylogenetic analysis of genomic or metagenomic sequence data";
}

sub usage_desc { "all %o <sequence file> [read 2 sequence file]" }

sub options {
	return (
		[ "force|f",      "Overwrites a previous Phylosift run with the same file name"],
		[ "custom=s",     "Reads a custom marker list from a file otherwise use all the markers from the markers directory"],
		[ "threads=i",    "Runs parallel portions using the specified number of processes (DEFAULT : 1)"],
		[ "extended",     "Uses the extended set of markers"],
		[ "paired",       "Looks for 2 input files (paired end sequencing) in FastQ format. Reversing the sequences for the second file listed and appending to the corresponding pair from the first file listed."],
		[ "continue",     "Enables the pipeline to continue to subsequent steps when not using the 'all' mode"],
		[ "isolate",      "Use this mode if you are running data from an isolate genome"],
		[ "simple",       "Creates a simple taxonomic summary of the output; no Krona output"],
		[ "besthit",      "When there are multiple hits to the same read, keeps only the best hit to that read"],
		[ "updated",      "Use the set of updated markers instead of stock markers"],
		[ "marker_url=s", "Phylosift will use markers available from the url provided"],
		[ "coverage=s",   "Provides a contig/scaffold coverage file to Phylosift"],
		[ "output=s",     "Specifies an output directory other than PStemp"],
		[ "stdin",        "Read sequence input on standard input"],
		[ "unique",       "Permit only a single hit between a marker and query sequence, discard any ambiguous hits"],
		[ "keep_search",  "Keeps the blastDir files (Default: Delete the blastDir files after every chunk)"],
		[ "chunk=i",      "Only run a set number of chunks"],
		[ "chunk_size=i", "Run so many sequences per chunk"],
	);
}

sub validate {
	my ($self, $opt, $args) = @_;

	# we need at least one argument beyond the options; die with that message
	# and the complete "usage" text describing switches, etc
	$self->usage_error("phylosift all requires exactly one or two file name arguments to run") unless @$args == 1 || @$args == 2;
}

sub load_opt {
	my %args = @_;
	my $opt = $args{opt};
	$Phylosift::Settings::file_dir = $opt->{output};
	$Phylosift::Settings::paired = $opt->{paired};
	$Phylosift::Settings::custom = $opt->{custom};
	$Phylosift::Settings::force = $opt->{force};
	$Phylosift::Settings::continue = $opt->{continue};
	$Phylosift::Settings::threads = $opt->{threads};
	$Phylosift::Settings::simple = $opt->{simple};
	$Phylosift::Settings::isolate = $opt->{isolate};
	$Phylosift::Settings::besthit = $opt->{besthit};
	$Phylosift::Settings::coverage = $opt->{coverage};
	$Phylosift::Settings::updated = $opt->{updated};
	$Phylosift::Settings::marker_url = $opt->{marker_url};
	$Phylosift::Settings::extended = $opt->{extended};
	$Phylosift::Settings::configuration = $opt->{config};
	$Phylosift::Settings::keep_search = $opt->{keep_search};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::unique = $opt->{unique};
	$Phylosift::Settings::start_chunk = $opt->{start_chunk};
	$Phylosift::Settings::stdin = $opt->{stdin};
	$Phylosift::Settings::chunks = $opt->{chunks};
	$Phylosift::Settings::my_debug = $opt->{debug};

	$Phylosift::Utilities::debuglevel = $Phylosift::Settings::my_debug || 0;
	
}

sub execute {
	my ($self, $opt, $args) = @_;
	load_opt{opt=>$opt};
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => $mode, file_1 => @$args[0], file_2 => @$args[1]);

	debug("FORCE: " . $Phylosift::Settings::force . "\n");
	debug("Continue : " . $Phylosift::Settings::continue . "\n");
	
	$ps->run( force=>$Phylosift::Settings::force, custom=>$Phylosift::Settings::custom, cont=>$Phylosift::Settings::continue );

}

1;