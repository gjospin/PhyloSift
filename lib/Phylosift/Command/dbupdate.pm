package Phylosift::Command::dbupdate;
use Phylosift -command;
use Phylosift::Settings;
use Phylosift::Phylosift;
use Phylosift::UpdateDB;
use Carp;
use Phylosift::Utilities qw(debug);

our $VERSION = "v1.0.1";
sub description {
	return "phylosift dbupdate - update the phylosift database with new genomic data";
}

sub abstract {
	return "update the phylosift database with new genomic data";
}

sub usage_desc { "dbupdate %o" }

sub options {
	return (
			 [ "repository=s",    "Path to repository for local copies of NCBI and EBI genome databases", ],
			 [ "destination=s",   "Path to destination for updated markers", ],
			 [ "knockouts=s",     "File containing a list of taxon IDs to exclude from the database" ],
			 [ "local-storage=s", "Path to local storage for results of scanning genomes for markers", ],
			 [ "base-markers=s",  "Path to base markers to add to updated set", ],
			 [ "skip-download",   "Skip downloading new data from external repositories" ],
			 [ "skip-scan",       "Skip scanning any new genomes for homologs to the marker database" ],
			 [ "keep-paralogs",   "Don't discard paralogs when building trees and fasta files" ],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;
}

sub set_ifdef {
	my $dest = $_[0];
	$$dest = $_[1] if defined $_[1];
}

=head2 check_required_opts

Checks to see if the 4 required options are defined before continuing with the update
They can be specified through the RC file or the command line

repository --repository
destination --destination
local_storage --local-storage
base_markers --base-markers

=cut

sub check_required_opts {
	my $error_message = "";
	$error_message .= "--repository needs to be specified and cannot be empty\n"
	  unless defined($Phylosift::Settings::repository) || $Phylosift::Settings::repository ne "";
	$error_message .= "--destination needs to be specified and cannot be empty\n"
	  unless defined($Phylosift::Settings::destination) || $Phylosift::Settings::repository ne "";
	$error_message .= "--local-storage needs to be specified and cannot be empty\n"
	  unless defined($Phylosift::Settings::local_storage) || $Phylosift::Settings::local_storage ne "";
	$error_message .= "--base-markers needs to be specified and cannot be empty\n"
	  unless defined($Phylosift::Settings::base_markers) || $Phylosift::Settings::base_markers ne "";
	croak($error_message) if $error_message ne "";
}

sub load_opt {
	my %args = @_;
	my $opt  = $args{opt};
	$Phylosift::Settings::configuration        = $opt->{config};
	$Phylosift::Settings::disable_update_check = $opt->{disable_updates};
	$Phylosift::Settings::my_debug             = $opt->{debug};
	$Phylosift::Settings::keep_paralogs        = $opt->{keep_paralogs};
	$Phylosift::Utilities::debuglevel          = $Phylosift::Settings::my_debug || 0;
}

sub execute {
	my ( $self, $opt, $args ) = @_;
	load_opt( opt => $opt );
	Phylosift::Command::sanity_check();
	my $ps = new Phylosift::Phylosift();
	set_ifdef( \$Phylosift::Settings::repository,      $opt->{repository} );
	set_ifdef( \$Phylosift::Settings::destination,     $opt->{destination} );
	set_ifdef( \$Phylosift::Settings::taxon_knockouts, $opt->{knockouts} );
	my $taxon_knockouts = $opt->{knockouts};

	#my $local = "/state/partition1/koadman/phylosift_extended_dbupdate/";
	set_ifdef( \$Phylosift::Settings::local_storage, $opt->{local_storage} );
	set_ifdef( \$Phylosift::Settings::base_markers,  $opt->{base_markers} );
	Phylosift::UpdateDB::set_default_values( post => 1 );
	check_required_opts();
	Phylosift::Utilities::program_checks();
	Phylosift::Utilities::data_checks( self => $ps );
	my $ebi_repository           = $Phylosift::Settings::repository."/ebi";
	my $ncbi_draft_repository    = $Phylosift::Settings::repository."/ncbi_draft";
	my $ncbi_finished_repository = $Phylosift::Settings::repository."/ncbi_finished";
	my $ncbi_wgs_repository      = $Phylosift::Settings::repository."/ncbi_wgs";
	my $local_repository         = $Phylosift::Settings::repository."/local";
	my $result_repository        = $Phylosift::Settings::destination."/processed";
	my $marker_dir               = $Phylosift::Settings::destination."/markers";
	my $manual_download          = $Phylosift::Settings::repository."/manual_download";
	my $newObject                = new Phylosift::Phylosift();
	my @new_genomes              = ();

	$newObject->{"updated"}              = 0;    # don't default to updated markers
	$newObject->{"disable_update_check"} = 1;    # don't update markers!
	Phylosift::Phylosift::read_phylosift_config( self => $newObject );
	Phylosift::Utilities::data_checks( self => $newObject );

	# force use of the new NCBI data
	debug "Updating NCBI tree and taxon map...";
	$Phylosift::Utilities::ncbi_dir = "$repository/ncbi/";
	Phylosift::UpdateDB::update_ncbi_taxonomy( repository => $Phylosift::Settings::destination );
	debug "Starting download\n";
	if ( !defined( $opt->{skip_download} ) ) {
		Phylosift::UpdateDB::get_ebi_genomes( directory => $ebi_repository );
		Phylosift::UpdateDB::get_ncbi_draft_genomes( directory => $ncbi_draft_repository );
		Phylosift::UpdateDB::get_ncbi_finished_genomes( directory => $ncbi_finished_repository );
		Phylosift::UpdateDB::get_ncbi_wgs_genomes( directory => $ncbi_wgs_repository );
	}
	if ( !defined( $opt->{skip_scan} ) ) {
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $ncbi_finished_repository,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $ncbi_draft_repository,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $ncbi_wgs_repository,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $ebi_repository,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $local_repository,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);
		Phylosift::UpdateDB::find_new_genomes(
											   genome_directory  => $manual_download,
											   results_directory => $Phylosift::Settings::local_storage,
											   files             => \@new_genomes
		);

		Phylosift::UpdateDB::qsub_updates( local_directory => $Phylosift::Settings::local_storage, files => \@new_genomes );
	}
	Phylosift::UpdateDB::collate_markers(
										  local_directory => $Phylosift::Settings::local_storage,
										  marker_dir      => $marker_dir,
										  taxon_knockouts => $Phylosift::Settings::taxon_knockouts
	);

	#Phylosift::UpdateDB::join_trees( marker_dir => $marker_dir );
	#	Phylosift::UpdateDB::package_markers( marker_directory => $marker_dir, base_marker_directory => $Phylosift::Settings::base_markers);

}

1;
