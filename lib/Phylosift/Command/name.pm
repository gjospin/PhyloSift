package Phylosift::Command::name;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::Summarize;
use Carp;
use Phylosift::Utilities qw(debug);

sub description {
	return "phylosift name <filename>- Renames the files in the available directories";
}

sub abstract {
	return "Renames the sequence IDs for all the available directories to the original IDs";
}

sub usage_desc { "name" }


sub options {
	my @opts = ();
	return @opts;
}

sub validate {
	my ($self, $opt, $args) = @_;
	$self->usage_error("phylosift name requires at least the filename as an argument") if @$args ==0 ;
	$self->usage_error("phylosift name can only have the filename as an argument") if @$args > 1;
}

sub execute {
	my ($self, $opt, $args) = @_;
	Phylosift::Command::all::load_opt(opt=>$opt);
	$Phylosift::Settings::keep_search = 1;
	Phylosift::Command::sanity_check();

	my $ps = new Phylosift::Phylosift();
	$ps = $ps->initialize( mode => "name", file_1 => @$args[0]);


	#gather the number of chunks to rename
	my @lookup_files = glob($ps->{"blastDir"} . "/lookup_ID.*.tbl");
	warn "Cannot rename sequences, no lookup file was found\n" if(@lookup_files ==0);
	foreach my $file (@lookup_files){
		$file =~ m/lookup_ID\.(\d+)\.tbl/;
		my $chunk = $1;
		Phylosift::Summarize::rename_sequences(self=>$ps, chunk=>$chunk) if (@lookup_files > 0);
	}
}

1;