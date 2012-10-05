package Phylosift::Command::index;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Carp;
use Phylosift::Utilities qw(debug);

sub description {
	return "phylosift index - index a phylosift database after changes have been made";
}

sub abstract {
	return "index a phylosift database after changes have been made";
}

sub usage_desc { "index" }


sub options {
	my @opts = ();
	return @opts;
}

sub validate {
	my ($self, $opt, $args) = @_;
	
	$self->usage_error("phylosift index does not accept any arguments") if @$args > 0;
}

sub execute {
	my ($self, $opt, $args) = @_;
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	Phylosift::Utilities::program_checks();
	Phylosift::Utilities::data_checks( self => $ps );

	@markers = Phylosift::Utilities::gather_markers(
		self              => $self,
		path              => $Phylosift::Settings::marker_dir,
		force_gather      => 1,
		allow_missing_hmm => 1
	);
	Phylosift::Utilities::index_marker_db(self => $self, markers => \@markers, path => $Phylosift::Settings::marker_dir);
	if (-d $Phylosift::Settings::markers_extended_dir && $Phylosift::Settings::extended ){
		my @extended_markers = Phylosift::Utilities::gather_markers(self => $self,path => $Phylosift::Settings::markers_extended_dir,
			force_gather => 1, allow_missing_hmm => 1 );
	
		Phylosift::Utilities::index_marker_db(
			self    => $self,
			markers => \@extended_markers,
			path    => $Phylosift::Settings::markers_extended_dir
		  );
	}
}

1;