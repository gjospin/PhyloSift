package Phylosift::Command::align;
use Phylosift::Command::all;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";

sub description {
	return "phylosift align - align homologous sequences identified by 'search'";
}

sub abstract {
	return "align homologous sequences identified by 'search'";
}

sub usage_desc { "align %o <sequence file> [pair sequence file]" }

sub align_opts {
	my %opts = (
				 besthit => [ "besthit", "When there are multiple hits to the same read, keeps only the best hit to that read",      { default => 0 } ],
				 stdin   => [ "stdin",   "Read sequence input on standard input" ],
				 unique  => [ "unique",  "Permit only a single hit between a marker and query sequence, discard any ambiguous hits", { default => 0 } ],
				 isolate => [ "isolate", "Use this mode if you are running data from an isolate genome",                             { default => 0 } ],
	);
	return %opts;
}

sub options {
	my %opts = align_opts();
	%opts = ( Phylosift::Command::all::all_opts(), %opts );
	return values(%opts);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
	Phylosift::Command::all::validate(@_);
	Phylosift::Command::all::validate_subcommand( @_, mode => "align" );
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	Phylosift::Command::all::load_opt( opt => $opt );
	$Phylosift::Settings::keep_search = 1;
	Phylosift::Command::sanity_check();
	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => "align", file_1 => @$args[0], file_2 => @$args[1] );
	$ps->{"ARGV"} = \@ARGV;
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::start_chunk, value => 1 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::chunks,      value => 1 );
	$ps->run( force => $Phylosift::Settings::force, custom => $Phylosift::Settings::custom, cont => $Phylosift::Settings::continue );
}

1;
