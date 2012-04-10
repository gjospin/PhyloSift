package Phylosift::Simulations;
use warnings;
use strict;
use Carp;
use Phylosift::Utilities;
use Phylosift::MarkerBuild;
use File::Basename;

=head1 SUBROUTINES/METHODS

=Head2 Simulations module

Generates Reads files simulated from randomly picked genomes from the concat.updated.pruned.tre in the markers directory.
Using only the DEFAULT marker package.
Builds a simulated marker package. This package is the default from which the picked genomes have been removed

Picks to be performed : Random
						Top PD

Need as input : Genome directory from where to look for the genomes picked
				Number of genomes to pick for each picks (Default = 10)
				Number of reads to generate per file (Default = 100,000)
				

=cut

#read simulations parameters
my $params_ill_fa = "-read_dist 105 -insert_dist 400 normal 50 -md poly4 3e-3 3.3e-8 -mr 95 5 ";
my $params_ill_fq = "-fq 1 -ql 30 10 " . $params_ill_fa;
my $params_454    = "-read_dist 100 normal 10 -homopolymer_dist balzer ";

=head2 prep_simulation

Runs the various steps necessary to running simulations
 - Pick representative genomes
 - Generates the input files (Illumina Fasta, 454 Fasta, paired ends FastQ in 2 files, paired ends FastQ interleaved)
 - Compresses the data using gzip
=cut

sub prep_simulation {
	my %args        = @_;
	my $self        = $args{self} || miss("PS_object");
	my $pick_number = $args{pick} || 10;
	my $read_number = $args{reads} || 100000;
	my $genomes_dir = $args{genomes_dir} || miss("Genome Directory");
	my $marker_dir  = Phylosift::Utilities::get_data_path( data_name => "markers", data_path => $Phylosift::Settings::marker_path );
	my $rep_file = Phylosift::MarkerBuild::get_representatives_from_tree(
																		  tree             => $marker_dir . "/concat.updated.pruned.tre",
																		  target_directory => $self->{"fileDir"},
																		  cutoff           => 0.01
	);
	my @genome_ids = get_genome_ids_from_pda( self => $self, rep_file => $rep_file, gene_map => "$marker_dir/gene_ids.aa.txt" );
	my ( $top_knockouts, $random_knockouts ) = knockout_genomes( self => $self, list => \@genome_ids, pick_number => $pick_number );
	debug "TOP $top_knockouts\n";
	debug "RAND $random_knockouts\n";
	my $top_ko_gen_file  = $self->{"fileDir"} . "/top_knockout.genomes";
	my $rand_ko_gen_file = $self->{"fileDir"} . "/random_knockout.genomes";
	my @gen_files        = ( $top_ko_gen_file, $rand_ko_gen_file );
	gather_genomes( self => $self, genome_list => $top_knockouts,    genomes => $genomes_dir, target => $top_ko_gen_file );
	gather_genomes( self => $self, genome_list => $random_knockouts, genomes => $genomes_dir, target => $rand_ko_gen_file );
	debug "Finished gathering the genomes\n";
	my ( $core, $path, $ext ) = fileparse( $top_ko_gen_file, qr/\.[^.]*$/ );

	foreach my $gen_file (@gen_files) {
		next unless -e $gen_file;
		#generate simulated reads for each knockout set
		simulate_reads(
						self    => $self,
						input   => $gen_file,
						reads   => $read_number * 2,
						params  => $params_ill_fq,
						outname => $core . "_ill_fastq",
						outdir  => $path
		);
		simulate_reads(
						self    => $self,
						input   => $gen_file,
						reads   => $read_number,
						params  => $params_ill_fa,
						outname => $core . "_ill_fasta",
						outdir  => $path
		);
		simulate_reads(
						self    => $self,
						input   => $gen_file,
						reads   => $read_number,
						params  => $params_454,
						outname => $core . "_454_fasta",
						outdir  => $path
		);
	}
}

=head2 gather_genomes

Gathers the genomes from a list and prints them to a file
Input needs a genome list, a genome directory and a file destination

=cut

sub gather_genomes {
	my %args        = @_;
	my $self        = $args{self} || miss("PS_object");
	my $genome_list = $args{genome_list} || miss("Genomes list");
	my $genomes_dir = $args{genomes} || miss("Genomes Directory");
	my $target_file = $args{target} || miss("Target File");
	my $FH          = Phylosift::Utilities::ps_open($genome_list);
	while (<$FH>) {
		chomp($_);
		my $grep_cmd = "ls $genomes_dir | grep \'.$_.fasta\'";
		debug "GREPPING : " . $grep_cmd . "\n";
		my $grep = `ls $genomes_dir | grep '.$_.fasta'`;
		if ( $grep eq "" ) {
			warn "Warning : $_ was not found\n";
		} else {

			#debug "GOT $grep";
			my @file_names = split( /\n/, $grep );
			my $file = $file_names[0];
			$file = get_largest_file( file_list => \@file_names, directory => $genomes_dir ) if scalar(@file_names) > 1;

			#debug "Using $file\n";
			`cat $genomes_dir/$file >> $target_file`;
		}
	}
}

=head2 get_largest_file

=cut

sub get_largest_file {
	my %args          = @_;
	my $file_list_ref = $args{file_list} || miss("File list reference");
	my $dir           = $args{directory} || miss("Directory");
	my $top_file      = "";
	my $top_size      = 0;
	foreach my $file ( @{$file_list_ref} ) {
		my $size = -s $dir . "/" . $file;
		if ( $size > $top_size ) {
			$top_file = $file;
			$top_size = $size;
		}
	}
	return $top_file;
}

=head2 simulate_reads

Simulates #### reads using Grinder from a genome list (Default : 100,000 reads generated)
Outputs the reads files in the PS_object's fileDir

Illumina Fasta, 454 Fasta, paired ends FastQ in 2 files, paired ends FastQ interleaved

This function required Grinder to be in the user's path

=cut

sub simulate_reads {
	my %args               = @_;
	my $self               = $args{self} || miss("PS_object");
	my $input_genomes_file = $args{input} || miss("Genomes list");
	my $read_number        = $args{reads} || 100000;
	my $out_directory      = $args{outdir} || miss("Output Directory");
	my $out_file           = $args{outname} || miss("Output file name");
	my $params             = $args{params} || miss("Simulate reads parameters");

	#Illumina paired ends  Generates
	debug "Simulating reads\n";
	my $simulation_cmd = "grinder " . $params . "-total_reads " . $read_number . " -bn $out_file -od $out_directory" . " -reference_file " . $input_genomes_file;
	debug "RUNNING : $simulation_cmd\n";
	`$simulation_cmd`;
}

=head2 get_genome_ids_from_pda

	Reads in a representatives file (output from pda) and compares them to the genes_ids.aa.txt in the marker directory to get genome_IDs
	Returns an array of genome_IDs.

=cut

sub get_genome_ids_from_pda {
	my %args     = @_;
	my $self     = $args{self} || miss("PS object");
	my $pda_file = $args{rep_file} || miss("pda_file");
	my $gene_map = $args{gene_map} || miss("Gene map");
	my %map      = ();
	my $FH       = Phylosift::Utilities::ps_open("$gene_map");
	while (<$FH>) {
		chomp($_);
		$_ =~ m/^(\S+)\s+(\S+)\s+(\S+)$/;
		$map{$3} = $2;
	}
	close($FH);

	#read in the available gene_oids
	#reading the pda file to get the representative IDs
	my $REPSIN = Phylosift::Utilities::ps_open("$pda_file");
	my @taxa   = ();
	while (<$REPSIN>) {
		chomp($_);
		next unless $_ =~ m/^(\d+)$/;
		if ( exists $map{$1} ) {
			push( @taxa, $map{$1} );
		} else {
			warn "Taxa not found\n";
		}
	}
	close($REPSIN);
	return @taxa;
}

=head2 knockout_genomes

	Prints 2 files in the fileDir for the PS object
		-top.ko contains a list of the knocked out taxa that have highest PD
		-random.ko contains a list of randomly knocked out taxa

=cut

sub knockout_genomes() {
	my %args        = @_;
	my $self        = $args{self} || miss("PS object");
	my $list_ref    = $args{list} || miss("Genome list reference");
	my $num_picked  = $args{pick_number} || 10;
	my @list        = @{$list_ref};
	my $list_length = scalar(@list);
	my $top_name    = $self->{"fileDir"} . "/top.ko";
	my $random_name = $self->{"fileDir"} . "/random.ko";
	my $TOP         = Phylosift::Utilities::ps_open(">$top_name");
	my $RAND        = Phylosift::Utilities::ps_open(">$random_name");
	foreach ( my $i = 0 ; $i < $num_picked ; $i++ ) {
		print $TOP $list[$i] . "\n";
		print $RAND $list[ int( rand($list_length) ) ] . "\n";
	}
	close($TOP);
	close($RAND);
	return ( $top_name, $random_name );
}
1;    # End of Phylosift::Simulations.pm
