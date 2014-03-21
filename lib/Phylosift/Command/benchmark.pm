package Phylosift::Command::benchmark;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::Simulations;
use File::Basename;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";

sub description {
	return "phylosift benchmark - measure taxonomic prediction accuracy on a simulated dataset";
}

sub abstract {
	return "measure taxonomic prediction accuracy on a simulated dataset";
}

sub usage_desc { "benchmark %o" }

sub options {
	return ( [ "output=s", "Path to write benchmark results", { default => "./" } ],
			 [ "summary-file=s", "Taxonomic summary input file", ],
			 [ "curve_path=s", "Path to write precision-recall curve data", { default => "./" } ],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;

	$self->usage_error("benchmark requires a phylosift output directory") unless @$args > 0;
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
	$ps->{"fileDir"} = "PS_temp/".@$args[0];    # this is the default location after running `phylosift all`
	my $summary = "PS_temp/".basename( @$args[0] )."/sequence_taxa.txt";
	$summary = $opt->{summary_file} if defined( $opt->{summary_file} );

	Phylosift::Benchmark::run_benchmark(
										 self         => $ps,
										 reads_file   => @$args[0],
										 output_path  => $opt->{output},
										 summary_file => $summary,
										 pr_curve     => $opt->{curve_path}
	);
}

1;
