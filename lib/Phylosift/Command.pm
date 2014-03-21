package Phylosift::Command;
use App::Cmd::Setup -command;
use Phylosift;
use POSIX;
use FindBin qw($Bin);
BEGIN { unshift( @INC, "$FindBin::Bin/../legacy/" ) if $] < 5.01; }

our $VERSION = "v1.0.1";

sub description {
	return "Analyze sequence data";
}

sub opt_spec {
	my ( $class, $app ) = @_;
	return (
		[ 'help' => "This usage screen" ],
		[ "debug",           "Print debugging messages" ],
		[ "disable_updates", "Disables automated check and download of Phylosift databases" ],
		[ "config=s",        "Provides a custom configuration file to Phylosift" ],
		$class->options($app),
	);
}

sub config {
	my $app = shift;
	Phylosift::read_phylosift_config();
}

sub sanity_check {
	unless ( (POSIX::uname)[4] =~ /64/ || $^O =~ /arwin/ ) {
		print STDERR (POSIX::uname)[4]."\n";
		die "Sorry, PhyloSift requires a 64-bit OS to run.\n";
	}
}

sub validate_args {
	my ( $self, $opt, $args ) = @_;
	if ( $opt->{help} ) {
		my ($command) = $self->command_names;
		$self->app->execute_command( $self->app->prepare_command( "help", $command ) );
		exit;
	}
	$self->validate( $opt, $args );
}

1;
