package Phylosift::Simulations;
use warnings;
use strict;
use Carp;
use Phylosift::Utilities;
use Phylosift::MarkerBuild;
use File::Basename;

our $VERSION = "v1.0.1";

=head1 SUBROUTINES/METHODS

=head2 Simulations module

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
my $params_ill_fa = "-read_dist 105 -insert_dist 300 normal 50 -md poly4 3e-3 3.3e-8 -mr 95 5 ";
my $params_ill_fq = "-fq 1 -ql 30 10 ".$params_ill_fa;
my $params_454    = "-read_dist 400 normal 40 -homopolymer_dist balzer ";

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

	# this finds the genomes which add the most PD to the tree
	my $rep_file = Phylosift::MarkerBuild::get_representatives_from_tree(
																		  tree             => $marker_dir."/concat.updated.tre",
																		  target_directory => $self->{"fileDir"},
																		  cutoff           => 0.01
	);

	# read a list of the available genomes for taxa
	debug "getting genome list\n";
	my $taxon_genome_files = find_taxon_genome_files( self => $self, genomes => $genomes_dir );

	my $new_genomes = pick_new_genomes( genomes => $genomes_dir, wgs_dir => "/share/eisen-d2/amphora2/ncbi_wgs/ftp.ncbi.nih.gov/genbank/wgs/" );

	debug "found ".scalar( keys(%$taxon_genome_files) )." genomes\n";

	# this reads in a list of genomes ordered by PD contribution
	my @genome_ids = get_genome_ids_from_pda( self => $self, rep_file => $rep_file, gene_map => "$marker_dir/gene_ids.aa.txt" );
	debug "Got ".scalar(@genome_ids)." genomes from PDA\n";
	debug "picking simulation genomes\n";
	my ($knockouts) =
	  knockout_genomes(
						self               => $self,
						list               => \@genome_ids,
						new_genomes        => $new_genomes,
						pick_number        => $pick_number,
						taxon_genome_files => $taxon_genome_files
	  );
	debug "knockouts $knockouts\n";
	my $ko_genome_file = $self->{"fileDir"}."/knockout.genomes";
	my %gen_files      = ( top => $ko_genome_file );
	my %gen_lists      = ( top => $knockouts );

	# now create a fasta file with the target genomes
	gather_genomes( self => $self, genome_list => $knockouts, genomes => $genomes_dir, target => $ko_genome_file );
	debug "Finished gathering the genomes\n";
	my ( $core, $path, $ext ) = fileparse( $ko_genome_file, qr/\.[^.]*$/ );

	foreach my $type ( keys(%gen_files) ) {
		next unless -e $gen_files{$type};

		#generate simulated reads for each knockout set
		simulate_reads(
						self              => $self,
						input             => $gen_files{$type},
						reads             => $read_number * 2,
						params            => $params_ill_fq." -am uniform 1",
						outname           => $core."unif_ill_fastq",
						outdir            => $path,
						distribution_file => $gen_lists{$type}
		);
		simulate_reads(
						self              => $self,
						input             => $gen_files{$type},
						reads             => $read_number * 2,
						params            => $params_ill_fq." -am exponential",
						outname           => $core."exp_ill_fastq",
						outdir            => $path,
						distribution_file => $gen_lists{$type}
		);
		simulate_reads(
						self              => $self,
						input             => $gen_files{$type},
						reads             => $read_number,
						params            => $params_454." -am uniform 1",
						outname           => $core."unif_454_fasta",
						outdir            => $path,
						distribution_file => $gen_lists{$type}
		);
		simulate_reads(
						self              => $self,
						input             => $gen_files{$type},
						reads             => $read_number,
						params            => $params_454." -am exponential",
						outname           => $core."exp_454_fasta",
						outdir            => $path,
						distribution_file => $gen_lists{$type}
		);
	}
}

sub find_taxon_genome_files {
	my %args        = @_;
	my $self        = $args{self} || miss("PS_object");
	my $genomes_dir = $args{genomes} || miss("Genomes Directory");
	my $LSSER       = ps_open("ls -1 $genomes_dir |");
	my %taxon_file_map;
	while ( my $line = <$LSSER> ) {
		chomp $line;
		my $taxon = $1 if $line =~ /\.(\d+)\.fasta$/;
		$taxon_file_map{$taxon} = $line if defined($taxon);
	}
	return \%taxon_file_map;
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
		my @line = split( /\t/, $_ );
		my $grep_cmd = "ls $genomes_dir | grep \'.$line[0].fasta\'";

		#debug "GREPPING : " . $grep_cmd . "\n";
		my $grep = `ls $genomes_dir | grep '.$line[0].fasta'`;
		if ( $grep eq "" ) {
			warn "Warning : $line[0] was not found\n";
		} else {

			#debug "GOT $grep";
			my @file_names = split( /\n/, $grep );
			my $file = $file_names[0];
			$file = get_largest_file( file_list => \@file_names, directory => $genomes_dir ) if scalar(@file_names) > 1;

			#debug "Using $file\n";
			my $OUT = Phylosift::Utilities::ps_open(">>$target_file");
			my $IN  = Phylosift::Utilities::ps_open("$genomes_dir/$file");
			while (<$IN>) {
				chomp($_);
				$_ =~ s/^>(\S+)/>$line[0] $1/g;
				print $OUT $_."\n";
			}
			close($OUT);
			close($IN);
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
		my $size = -s $dir."/".$file;
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
	my $distrib            = $args{distribution_file} || miss("Distribution file");

	#Illumina paired ends  Generates
	debug "Simulating reads\n";
	my $simulation_cmd = "grinder ".$params." -total_reads ".$read_number." -bn $out_file -od $out_directory"." -reference_file ".$input_genomes_file;
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

# /share/eisen-d2/amphora2/ncbi_wgs/ftp.ncbi.nih.gov/genbank/wgs/

sub pick_new_genomes {
	my %args        = @_;
	my $genomes_dir = $args{genomes} || miss("Genomes Directory");
	my $wgs_dir     = $args{wgs_dir} || miss("WGS Directory");
	my %new_genomes = ();
	my $LSSER       = ps_open("ls -1 $genomes_dir |");
	my %taxon_file_map;
	while ( my $line = <$LSSER> ) {
		chomp $line;
		my ( $wgs_id, $taxon_id ) = ( $1, $2 ) if $line =~ /(....)\.(\d+)\.fasta$/;
		my $gbff = "$wgs_dir/wgs.$wgs_id.1.gbff";
		$gbff = "$wgs_dir/wgs.$wgs_id.gbff" unless -e $gbff;
		next unless -e $gbff;
		my $WGS_IN = ps_open($gbff);
		my $fline  = <$WGS_IN>;
		my @dat    = split( /\s+/, $fline );
		my @date   = split( /-/, $dat[7] );
		next unless $date[2] >= 2012;
		next unless $date[1] eq "APR" || $date[1] eq "MAY" || $date[1] eq "JUN" || $date[1] eq "MAR";
		$new_genomes{$taxon_id} = 1;
	}
	print STDERR "Found ".scalar( keys(%new_genomes) )." new genomes\n";
	return \%new_genomes;
}

=head2 knockout_genomes

	Prints 2 files in the fileDir for the PS object
		-top.ko contains a list of the knocked out taxa that have highest PD
		-random.ko contains a list of randomly knocked out taxa

=cut

sub knockout_genomes() {
	my %args               = @_;
	my $self               = $args{self} || miss("PS object");
	my $list_ref           = $args{list} || miss("Genome list reference");
	my $num_picked         = $args{pick_number} || 10;
	my $taxon_genome_files = $args{taxon_genome_files};
	my $new_genomes        = $args{new_genomes};
	my @list               = @{$list_ref};
	my $list_length        = scalar(@list);
	my $ko_name            = $self->{"fileDir"}."/grinder_knockouts.txt";
	my $full_ko_name       = $self->{"fileDir"}."/knockouts.txt";
	my $KO                 = Phylosift::Utilities::ps_open(">$ko_name");
	my $seed_percent       = 100;
	my @top_array          = @list[ 0 .. $num_picked - 1 ];
	my %ko_taxa;
	my $abund_sum = 0;
	my $nummer    = $num_picked;
	my $use_maxpd = 0;

	$num_picked *= 2 unless $use_maxpd;    # gawd this is ugly.

	#	for ( my $i = 0 ; $i < $num_picked/2 ; ) {
	#		# pick one from the maxpd set
	#		my $rand_index = int( rand( scalar(@top_array) - 1 ) );
	#		my $rand_taxon  = $top_array[$rand_index];
	#		$top_array[$rand_index] = $list[$nummer++];
	#			debug "picked $rand_taxon, i $i\n";
	#		next unless defined($taxon_genome_files->{$rand_taxon}); # don't add this one unless its available
	#		$ko_taxa{$rand_taxon}=rand(1000);	# assign a uniformly random abundance
	#		$abund_sum += $ko_taxa{$rand_taxon};
	#		splice( @top_array, $rand_index, 1 );	# remove this one so it doesn't get resampled
	#		$i++;
	#	}
	#	for ( my $i = $num_picked/2 ; $i < $num_picked ; ) {
	my @newgenomes = keys(%$new_genomes);
	for ( my $i = 0; $i < $num_picked; ) {

		# pick one uniformly at random
		#		my $rand_taxon = $list[ int( rand($list_length) ) ];
		my $rand_taxon = $newgenomes[ int( rand(@newgenomes) ) ];
		debug "picked $rand_taxon, i $i\n";
		croak("ran out of genomes, relax your constraints!") if $nummer == @list;
		print "$rand_taxon not found on disk " unless defined( $taxon_genome_files->{$rand_taxon} );
		print "$rand_taxon not new "           unless defined( $new_genomes->{$rand_taxon} );
		next                                   unless defined( $taxon_genome_files->{$rand_taxon} );    # is it available on disk?
		next                                   unless defined( $new_genomes->{$rand_taxon} );           # is it new enough?
		next if defined( $ko_taxa{$rand_taxon} );                                                       # do we already have this one?
		$ko_taxa{$rand_taxon} = rand(1000);
		$abund_sum += $ko_taxa{$rand_taxon};
		$i++;
	}
	foreach my $taxon ( keys(%ko_taxa) ) {
		print $KO "$taxon\n";                                                                           #.($ko_taxa{$taxon}/$abund_sum)."\n";
	}
	close($KO);
	my @kos      = keys(%ko_taxa);
	my $extra_ko = find_neighborhoods( self => $self, taxa => \@kos, swath => 0.15, concat_tree => "$Phylosift::Utilities::marker_dir/concat.updated.tre" );
	my $FULLKO   = Phylosift::Utilities::ps_open(">$full_ko_name");
	foreach my $taxon (@$extra_ko) {
		print $FULLKO "$taxon\n";
	}

	return $ko_name;
}

sub find_neighborhoods {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $concat_tree = $args{concat_tree} || miss("concat_tree");
	my $taxa        = $args{taxa} || miss("taxa");
	my $swath       = $args{swath} || miss("swath");

	debug "reading concat tree\n";
	my $tree = Bio::Phylo::IO->parse(
									  '-file'   => $concat_tree,
									  '-format' => 'newick',
	)->first;

	# find the taxon of interest
	debug "reading gene IDs\n";
	my ( $id_to_taxon, $marker_taxon_to_id ) = Phylosift::UpdateDB::read_gene_ids( file => "$Phylosift::Utilities::marker_dir/gene_ids.aa.txt" );
	my @target_nodes;
	my %del_nodes;
	foreach my $taxon (@$taxa) {
		$del_nodes{$taxon} = 1;
		my $gene_id = $marker_taxon_to_id->{"concat"}{$taxon}->[0];
		my $target_node;
		debug "Looking for neighbors of $gene_id\n";
		foreach my $node ( @{ $tree->get_entities } ) {
			if ( $node->get_name eq $gene_id ) {
				$target_node = $node;
				push( @target_nodes, $node );
				last;
			}
		}
		croak("Unable to find node for $gene_id, taxon $taxon") unless defined $target_node;

		# find all neighbors within distance X
		foreach my $node ( @{ $tree->get_entities } ) {

			# skip this one if it is not a leaf
			next if ( scalar( @{ $node->get_children() } ) > 1 );
			my $d1 = $node->calc_patristic_distance($target_node);
			next if ( $d1 > $swath );
			my $name = $node->get_name;
			$del_nodes{ $id_to_taxon->{$name} } = 1;
		}
	}

	# now calculate stats.
	my %valid;    # hash of nodes still in tree
	foreach my $node ( @{ $tree->get_entities } ) {

		# for each non-deleted leaf, walk up the tree and record it as valid
		next if ( scalar( @{ $node->get_children() } ) > 1 );
		next if defined( $del_nodes{ $id_to_taxon->{ $node->get_name } } );
		my $parent = $node;
		while ( $parent = $parent->get_parent() ) {
			$valid{$parent} = $parent;
		}
	}
	debug "Started with ".scalar( @{ $tree->get_entities } )." nodes\n";
	debug scalar( keys(%valid) )." valid nodes remaining\n";

	# each valid node should have zero or >= 2 valid children
	my %dvalid;
	foreach my $np ( keys(%valid) ) {
		my $v_child = 0;
		foreach my $child ( @{ $valid{$np}->get_children() } ) {
			$v_child++ if ( defined( $valid{$child} ) );
		}
		$dvalid{$np} = $valid{$np} if $v_child != 1;
	}
	debug scalar( keys(%dvalid) )." doubly valid nodes remaining\n";

	# now for each deleted node, find it's ancestral valid node & record some stats
	my %stats;
	my $SIMSTATS = ps_open( ">".$self->{"fileDir"}."/sim_stats.txt" );
	print $SIMSTATS
	  "# TaxonID\tParentID\tParentDistance\tParentToGrandParent\tParentOtherChild\tMinLeafDist\tNodes0.05\tNodes0.10\tNodes0.15\tNodes0.20\tNodes0.25\n";

	foreach my $node (@target_nodes) {
		my $parent = $node;
		my $p_len  = $node->get_branch_length();
		while ( defined( $parent->get_parent() ) ) {
			$parent = $parent->get_parent();
			last if $dvalid{$parent};
			$p_len += $parent->get_branch_length() if defined( $parent->get_branch_length() );
		}
		debug "dist to dvalid parent is $p_len\n";

		# get distance from ancestor to grandpere
		my $pp_len = defined( $parent->get_branch_length() ) ? $parent->get_branch_length() : 0;
		debug "pplen $pp_len\n";

		# get distance to other valid child
		my $c_dist = 0;
		my $c      = $parent;
		while ( defined($c) && ( $c == $parent || !defined( $dvalid{$c} ) ) ) {
			my $new_c;
			foreach my $child ( @{ $c->get_children() } ) {
				next unless defined( $valid{$c} );
				$c_dist += $c->get_branch_length() if defined( $c->get_branch_length() );
				$new_c = $child;
				last;
			}
			$c = $new_c;
		}

		# get distance to nearest valid leaf
		my $min_leaf_dist = -1;
		foreach my $lnode ( @{ $tree->get_entities } ) {
			next unless ( @{ $lnode->get_children() } < 2 );
			my $name = $lnode->get_name;
			next if defined( $del_nodes{ $id_to_taxon->{$name} } );
			my $ldist = $lnode->calc_patristic_distance($node);
			$min_leaf_dist = $ldist if ( $min_leaf_dist < 0 || $ldist < $min_leaf_dist );
		}

		# calculate node densities
		my %n_density;
		for ( my $ddist = 5; $ddist <= 50; $ddist += 5 ) {
			$n_density{$ddist} = 0;
			foreach my $lnode ( @{ $tree->get_entities } ) {
				next unless $dvalid{$lnode};
				next unless $lnode->calc_patristic_distance($node) < $ddist / 100;
				$n_density{$ddist}++;
			}
		}

		my $tid = $id_to_taxon->{ $node->get_name() };
		$stats{$node} = [ $tid, $parent, $p_len, $pp_len, $c_dist ];

		# seems like there's some problem with c_dist
		print $SIMSTATS join( "\t",
							  $tid,           $parent,        $p_len,         $pp_len,        $c_dist,        $min_leaf_dist, $n_density{5},  $n_density{10},
							  $n_density{15}, $n_density{20}, $n_density{25}, $n_density{30}, $n_density{35}, $n_density{40}, $n_density{45}, $n_density{50} )
		  ."\n";
	}
	my @dn = keys(%del_nodes);
	return \@dn;
}

1;    # End of Phylosift::Simulations.pm
