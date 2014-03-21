package Phylosift::Command::summarize;
use Phylosift::Command::all;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";

sub description {
	return "phylosift summarize - translate a collection of phylogenetic placements into a taxonomic summary";
}

sub abstract {
	return "translate a collection of phylogenetic placements into a taxonomic summary";
}

sub usage_desc { "summarize %o <sequence file> [pair sequence file]" }

sub summarize_opts {
	my %opts = ( simple => [ "simple", "Creates a simple taxonomic summary of the output; no Krona output" ], );
	return %opts;
}

sub options {
	my %opts = summarize_opts();
	%opts = ( Phylosift::Command::all::all_opts(), %opts );
	return values(%opts);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
	Phylosift::Command::all::validate(@_);
	Phylosift::Command::all::validate_subcommand( @_, mode => "summarize" );
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	Phylosift::Command::all::load_opt( opt => $opt );
	$Phylosift::Settings::keep_search = 1;
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => "summarize", file_1 => @$args[0], file_2 => @$args[1] );
	$ps->{"ARGV"} = \@ARGV;
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::start_chunk, value => 1 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::chunks,      value => 1 );
	$ps->run( force => $Phylosift::Settings::force, custom => $Phylosift::Settings::custom, cont => $Phylosift::Settings::continue );
}

1;
