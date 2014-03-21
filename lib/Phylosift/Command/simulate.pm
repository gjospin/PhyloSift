package Phylosift::Command::simulate;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::Simulations;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";
sub description {
	return qq{phylosift simulate - simulate sequencing from a metagenomic sample
Example: 
> phylosift.pl sim --genome_dir=genomes --genome_count=10 --read_count=100000
};
}

sub abstract {
	return "simulate sequencing from a metagenomic sample";
}

sub usage_desc { "simulate %o" }

sub options {
	return (
			 [ "genome_count=i", "The number of genomes to include in the simulated community", { default => 10 } ],
			 [ "read_count=i",   "Number of reads in the simulated community",                  { default => 100000 } ],
			 [ "genome_dir=s",   "Path to a genome repository created by phylosift updateDB" ],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
}

sub load_opt {
	my %args = @_;
	my $opt  = $args{opt};
	$Phylosift::Settings::configuration        = $opt->{config};
	$Phylosift::Settings::keep_search          = $opt->{keep_search};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::my_debug             = $opt->{debug};

	$Phylosift::Utilities::debuglevel = $Phylosift::Settings::my_debug || 0;
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	load_opt( opt => $opt );
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	Phylosift::Simulations::prep_simulation( self => $ps, pick => $opt->{genome_count}, reads => $opt->{read_count}, genomes_dir => $opt->{genome_dir} );
}

1;
