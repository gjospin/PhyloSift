package Phylosift::Command::search;
use Phylosift::Command::all;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::FastSearch;
use Carp;
use Phylosift::Utilities qw(debug);

sub description {
	return "phylosift search - search input sequence for homology to reference gene database";
}

sub abstract {
	return "search input sequence for homology to reference gene database";
}

sub usage_desc { "search %o <sequence file> [pair sequence file]" }

sub search_opts {
	my %opts = (
		besthit => [ "besthit",      "When there are multiple hits to the same read, keeps only the best hit to that read", {default => 0}],
		stdin =>   [ "stdin",        "Read sequence input on standard input"],
		unique =>  [ "unique",       "Permit only a single hit between a marker and query sequence, discard any ambiguous hits", {default => 0}],
		isolate => [ "isolate",      "Use this mode if you are running data from an isolate genome", {default => 0}],
	);
	return %opts;
}

sub options {
	my %opts = search_opts();	
	%opts = (Phylosift::Command::all::all_opts(), %opts);
	return values(%opts);
}

sub validate {
	my ($self, $opt, $args) = @_;
	
	$self->usage_error("phylosift search requires exactly one or two file name arguments to run") unless @$args == 1 || @$args == 2;
}

sub execute {
	my ($self, $opt, $args) = @_;
	Phylosift::Command::all::load_opt(opt=>$opt);
	$Phylosift::Settings::keep_search = 1;
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => "search", file_1 => @$args[0], file_2 => @$args[1]);
	$ps->{"ARGV"} = \@ARGV;
	$ps->run( force=>$Phylosift::Settings::force, custom=>$Phylosift::Settings::custom, cont=>$Phylosift::Settings::continue );
}

1;