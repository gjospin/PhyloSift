package Phylosift::Command::name;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::Summarize;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";

sub description {
	return "phylosift name <filename>- Swap the original sequence names back into output files";
}

sub abstract {
	return "Replaces phylosift's own sequence IDs with the original IDs found in the input file header";
}

sub usage_desc { "name <sequence file> [pair sequence file]" }

sub name_opts {
	my %opts = (
				 besthit => [ "besthit", "When there are multiple hits to the same read, keeps only the best hit to that read", { default => 0 } ],
				 isolate => [ "isolate", "Use this mode if you are running data from an isolate genome",                        { default => 0 } ],
	);
	return %opts;
}

sub options {
	my %opts = name_opts();
	%opts = ( Phylosift::Command::all::all_opts(), %opts );
	return values(%opts);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
	$self->usage_error("phylosift name requires at least the filename as an argument") if @$args == 0;
	$self->usage_error("phylosift name can only have the filename(s) as an argument")  if @$args > 2;
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	Phylosift::Command::all::load_opt( opt => $opt );
	$Phylosift::Settings::keep_search = 1;
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => "name", file_1 => @$args[0] );

	#gather the number of chunks to rename
	my @lookup_files = glob( $ps->{"blastDir"}."/lookup_ID.*.tbl" );
	warn "Cannot rename sequences, no lookup file was found\n" if ( @lookup_files == 0 );
	foreach my $file (@lookup_files) {
		$file =~ m/lookup_ID\.(\d+)\.tbl/;
		my $chunk = $1;
		Phylosift::Summarize::rename_sequences( self => $ps, chunk => $chunk ) if ( @lookup_files > 0 );
	}
}

1;
