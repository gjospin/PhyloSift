package Phylosift::MarkerBuild;
use Cwd;
use Carp;
use Phylosift::Utilities qw(:all);
use File::Basename;
use JSON;

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::MarkerBuild - build a seed marker package from an existing multiple sequence alignment

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS

Given a multiple sequence alignment, functions in this package can build a marker package.
Steps involved include masking the alignment, selecting representative sequences for database searches, building a phylogenetic tree, HMM, and pplacer package.

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 build_marker

=cut
sub build_marker {
	my %args        = @_;
	my $opt         = $args{opt};
	my $self        = $args{self} || miss("self");
	my $aln_file    = $args{alignment} || miss("alignment");
	my $reps_pd     = $args{reps_pd} || miss("reps_pd");
	my $tree_pd     = $args{tree_pd} || miss("tree_pd");
	my $force       = $args{force};                            #force is not a required option
	my $mapping     = $args{mapping};                          #not a required argument
	my $destination = $args{destination};
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $marker_dir = Phylosift::Utilities::get_data_path( data_name => "markers",
														  data_path => $Phylosift::Settings::marker_path );

	$destination = $marker_dir unless defined($destination);
	$target_dir = $destination."/$core";

	# check that taxit is available
	my $taxit = Phylosift::Utilities::get_program_path( prog_name => "taxit" );
	if ( $taxit eq "" ) {
		croak(
			"Error: you must install pplacer's taxtastic and its dependencies prior to building markers. See https://github.com/fhcrc/taxtastic for more information.\n"
		);
	}

	# check for other programs
	Phylosift::Utilities::program_checks();

	# load NCBI taxonomy and marker info
	Phylosift::Utilities::data_checks( self => $self );

	#$target_dir = $Phylosift::Settings::file_dir;
	if ( -e $target_dir && !$force ) {
		croak(
			"Marker already exists in $destination. Delete Marker and restart the marker build.\nUse -f to force an override of the previous marker.\nUsage:\n>phylosift build_marker -f aln_file cutoff\n"
		);
	} else {
		`rm -rf "$target_dir"` if $force;
		`mkdir "$target_dir"`;
	}
	my $fasta_file = "$target_dir/$core.fasta";
	$aln_file = check_sequence_integrity( input => $aln_file, output_dir => $target_dir );

	my $clean_aln = $aln_file;
	my $seq_count = -1;
	unless ( $opt->{update_only} ) {

		# this code path generates new alignments
		# no need to do this if we're just updating an existing marker with new sequences
		$seq_count = Phylosift::Utilities::unalign_sequences( aln => $aln_file, output_path => $fasta_file );
		my $masked_aln = "$target_dir/$core.masked";
		mask_aln( file => $aln_file, output => $masked_aln ) unless -e $masked_aln;
		my $hmm_file = "$target_dir/$core.hmm";
		generate_hmm( file_name => $aln_file, hmm_name => $hmm_file ) unless -e $hmm_file;
		Phylosift::Utilities::fasta2stockholm( fasta => $aln_file, output => "$target_dir/$core.stk" ) unless -e "$target_dir/$core.stk";
		my $stk_aln = "$target_dir/$core.stk";

		#may need to create an unaligned file for the sequences before aligning them
		my $new_alignment_file = hmmalign_to_model(
													hmm_profile         => $hmm_file,
													sequence_file       => $fasta_file,
													target_dir          => $target_dir,
													reference_alignment => $stk_aln,
													sequence_count      => $seq_count
		);
		$aln_file = $new_alignment_file;
	} elsif ( -e $opt->{unaligned} ) {

		# copy the unmasked sequences if they exist
		debug "Copying the unaligned fasta sequences\n";
		`cp $opt->{unaligned} $fasta_file`;
	}
	$clean_aln = "$target_dir/$core.clean";
	my %id_map = mask_and_clean_alignment( alignment_file => $aln_file, output_file => $clean_aln );
	debug( "ID_map is ".scalar( keys(%id_map) )." long\n" );
	my ( $fasttree_file, $tree_log_file ) = generate_fasttree( alignment_file => $clean_aln, target_directory => $target_dir )
	  unless -e "$target_dir/$core.tree";
	pd_prune_fasta( tre => $fasttree_file, distance => $tree_pd, fasta => $clean_aln, pruned_fasta => "$target_dir/$core.pruned.fasta" );

	( $fasttree_file, $tree_log_file ) = generate_fasttree( alignment_file => "$target_dir/$core.pruned.fasta", target_directory => $target_dir );
	`mv $target_dir/$core.pruned.fasta $clean_aln`;
	debug "CLEAN_ALN : $clean_aln";

	#mask_and_clean_alignment( alignment_file => "$target_dir/$core.pruned.fasta", output_file => $clean_aln );
	my $rep_file;
	if ( $seq_count > 10 || $seq_count < 0 ) {
		debug "Looking for representatives\n";

		#need to generate representatives using PDA
		$rep_file = get_representatives_from_tree(
												   tree             => $fasttree_file,
												   target_directory => $target_dir,
												   reps_pd          => $reps_pd
		) unless -e "$target_dir/$core.pda";

		# need to read the representatives picked by PDA and generate a representative fasta file
		my $rep_fasta = get_fasta_from_pda_representatives(
															pda_file        => $rep_file,
															target_dir      => $target_dir,
															fasta_reference => $fasta_file,
															id_map          => \%id_map,
															core            => $core,
		) if -e $fasta_file;
	} else {

		#use all the sequences for representatives
		`cp $fasta_file $target_dir/$core.rep` if -e $fasta_file;
	}

	create_taxon_table( target_dir => $target_dir, mapping => $mapping, id_map => \%id_map ) if defined($mapping);

	if ( defined $mapping ) {
		my $tmp_jplace    = $target_dir."/".$core.".tmpread.jplace";
		my $mangled       = $target_dir."/".$core.".tmpread.mangled";
		my $taxon_map     = $target_dir."/".$core.".taxonmap";
		my $id_taxon_map  = $target_dir."/".$core.".gene_map";
		my $ncbi_sub_tree = $target_dir."/".$core.".subtree";
		debug("Using the mapping stuff\n");
		debug(
			"taxit create -a \"Guillaume Jospin\" -d \"simple package for reconciliation only\" -l temp -f $clean_aln -t $fasttree_file -s $tree_log_file -P $target_dir/temp_ref"
		);

		#make dummy placement reference package
		`taxit create -a "Guillaume Jospin" -d "simple package for reconciliation only" -l temp -f $clean_aln -t $fasttree_file -s $tree_log_file -P $target_dir/temp_ref`
		  unless -e "target_dir/temp_ref";

		#make a dummy jplace file
		my $tmpread_file = create_temp_read_fasta( file => "$target_dir/$core", aln_file => $clean_aln ) unless -e "$target_dir/$core.tmpread.fasta";
		`cd "$target_dir";$Phylosift::Settings::pplacer -c temp_ref -p "$tmpread_file"` unless -e $tmp_jplace;
		tree_mangler( in => $tmp_jplace, out => $mangled );

		#create a file with a list of IDs
		my %taxon_hash = generate_id_to_taxonid_map(
													 self        => $self,
													 id_hash_ref => \%id_map,
													 map_file    => $mapping,
													 out_file    => $id_taxon_map,
													 marker      => $core
		);
		my @taxon_array = keys(%taxon_hash);
		debug "TAXON_ARRAY ".scalar(@taxon_array)."\n";
		make_ncbi_subtree( out_file => $ncbi_sub_tree, taxon_ids => \@taxon_array );
		debug "AFTER ncbi_subtree\n";
		my $rconc = "$Phylosift::Settings::readconciler $ncbi_sub_tree $mangled $id_taxon_map $taxon_map";
		debug $rconc."\n";
		system($rconc);
		debug "AFTER readconciler\n";

		`rm -rf "$tmpread_file" "$mangled" "$tmp_jplace" "$ncbi_sub_tree"  "$target_dir/temp_ref"`;
	}

	#use taxit to create a new reference package required for running PhyloSift
	#needed are : 1 alignment file, 1 representatives fasta file, 1 hmm profile, 1 tree file, 1 log tree file.
	my $taxdb_opts = "";    # add taxit-friendly taxon labels if available
	                        #$taxdb_opts = "-i $target_dir/seq_ids.csv -T $target_dir/taxa.csv" if -e "$target_dir/taxa.csv" && -s "$target_dir/taxa.csv" > 0;
	my $taxit_cmd =
	  "cd \"$target_dir\";taxit create -c -d \"Creating a reference package for PhyloSift for the $core marker\" -l \"$core\" -f \"$clean_aln\" -t \"$fasttree_file\" $taxdb_opts -s \"$tree_log_file\" -P \"$core\"";
	debug "Running $taxit_cmd\n";
	`$taxit_cmd`;

	`rm -f "$target_dir/$core.pda"` if -e "$target_dir/$core.pda";
	`rm -f "$target_dir/$core.pda"` if -e "$target_dir/$core.pruned.pda";
	`rm -f "$target_dir/$core.tree"`;
	`rm -f "$target_dir/$core.log"`;
	`rm -f "$target_dir/$core.aln"`;
	`rm -f "$target_dir/$core.fasta"`;
	`rm -f "$target_dir/$core.checked"` if -e "$target_dir/$core.checked";
	`rm -f "$clean_aln"` unless $opt->{update_only};
	`mv -f "$target_dir/$core"/* "$target_dir"`;
	`rm -rf "$target_dir/$core"`;
	`rm -rf "$Phylosift::Settings::file_dir"`;
}

sub leaf_only_readconciler {
	my %args         = @_;
	my $mapping_file = $args{mapping} || miss("Missing mapping file\n");
	my $jplace       = $args{jplace} || miss("jplace file \n");
	my $output       = $args{output} || miss("output file\n");
	my $id_map_ref   = $args{id_map_ref} || miss("id map hash\n");
	my $JPLACEINPUT  = ps_open($jplace);
	my @treedata     = <$JPLACEINPUT>;
	close $JPLACEINPUT;
	my %unique_id_map         = %{$id_map_ref};
	my $json_data             = decode_json( join( "", @treedata ) );
	my $tree_string           = $json_data->{tree};
	my %unique_to_original_id = ();

	foreach my $key ( keys %unique_id_map ) {

		#print $key."\t".$unique_id_map{$key}."\n";
		$unique_to_original_id{ $unique_id_map{$key} } = $key;
	}
	print "Using MAPPING file : $mapping_file\n";
	$tree_string =~ s/^\s+\"//g;
	$tree_string =~ s/:(.+?)(\{\d+?\})/$2:$1/g;
	$tree_string =~ s/\"\,$//g;

	my $OUTPUT = ps_open(">$output");

	my %map     = ();
	my $MAPPING = ps_open($mapping_file);
	while (<$MAPPING>) {
		chomp($_);
		$_ =~ m/^(\S+)\s+(\d+)$/;
		$map{$1} = $2;
	}

	my $tree = Bio::Phylo::IO->parse( '-string' => $tree_string, '-format' => 'newick' )->first;
	foreach my $node ( @{ $tree->get_entities } ) {
		my $name = $node->get_name;
		next unless $name =~ m/^(\d+)\{(\d+)\}$/;
		if ( exists $map{ $unique_to_original_id{$1} } ) {
			print $OUTPUT "$2\t$map{$unique_to_original_id{$1}}\n";
		}
	}
	close($MAPPING);
	close($OUTPUT);
}

sub create_taxon_table {
	my %args       = @_;
	my $mapping    = $args{mapping};
	my $id_map     = $args{id_map};
	my $target_dir = $args{target_dir};

	# create a taxon table for taxit
	my $TAXIN      = ps_open($mapping);
	my $TAXIDTABLE = ps_open(">$target_dir/tax_ids.txt");
	my $SEQIDS     = ps_open(">$target_dir/seq_ids.csv");
	print $SEQIDS "\"seqname\",\"tax_id\"\n";
	my %allanc;
	while ( my $line = <$TAXIN> ) {
		chomp $line;
		my ( $seq_name, $tid ) = split( /\t/, $line );
		my $seq_id = $id_map->{$seq_name};
		next unless $tid;
		my @tinfo = Phylosift::Summarize::get_taxon_info( taxon => $tid );   # lookup will check for any merging
		                                                                     # not sure what taxit doesn't like about the astrovirus but it won't handle them...
		if ( length( $tinfo[2] ) > 0 && $tinfo[2] ne "ERROR" && $tinfo[0] !~ /ASTROVIRUS/ ) {
			$allanc{ $tinfo[2] } = 1;
		}
		print $SEQIDS "$seq_id,$tid\n";
	}
	print $TAXIDTABLE join( "\n", keys(%allanc) );
	close $TAXIDTABLE;
	my $t_cmd = "taxit taxtable -d $target_dir/../ncbi_taxonomy.db -t $target_dir/tax_ids.txt -o $target_dir/taxa.csv";
	debug "Running $t_cmd\n";
	system($t_cmd);
}

=head2 check_sequence_integrity

Checks for characters in the input fasta file. Removes/changes characters accordingly. Exits PS if bad characters are found.

=cut

sub check_sequence_integrity {
	my %args       = @_;
	my $file_input = $args{input} || miss("Input file to check");
	my $output_dir = $args{output_dir} || miss("Output directory");
	debug "Checking integrity on File $file_input\n";
	my ( $core, $path, $ext ) = fileparse( $file_input, qr/\.[^.]*$/ );
	my $type_ref = Phylosift::Utilities::get_sequence_input_type( ps_open($file_input) );
	my %type     = %{$type_ref};
	croak("Input was detected as DNA.\n Terminating\n") if $type{seqtype} eq 'DNA';
	my $IN_HANDLE  = ps_open($file_input);
	my $OUT_HANDLE = ps_open(">$output_dir/$core.checked");

	while (<$IN_HANDLE>) {

		if ( $_ =~ m/^>/ ) {

			#do nothing
		} elsif ( $_ =~ m/\*/ ) {
			if ( $type{seqtype} eq 'protein' ) {
				$_ =~ s/\*/X/g;
			} else {
				##should never reach this anymore but keeping this in in case we allow DNA in the Future.
				croak("Input was detected as DNA and a bad character was found.\nTerminating.\n;");
			}
		} elsif ( $_ =~ m/U/ && $type{seqtype} eq 'DNA' ) {
			##should never reach this anymore but keeping this in in case we allow DNA in the Future.
			croak("Input was detected as DNA and a Uracil character was found.\nTerminating.\n;");
		}
		print $OUT_HANDLE $_;
	}
	close($IN_HANDLE);
	close($OUT_HANDLE);
	return ("$output_dir/$core.checked");
}

=head2 create_temp_read_fasta

reads a fasta file and writes 1 read/sequence to the output file specified

=cut

sub create_temp_read_fasta {
	my %args     = @_;
	my $file     = $args{file} || miss("file");
	my $aln_file = $args{aln_file} || miss("aln_file");
	my $TMPREAD  = ps_open(">$file.tmpread.fasta");
	debug "creating tmpread $file.tmpread.fasta\n";

	print $TMPREAD ">blahblahblah\n";
	my $ALNIN       = ps_open("$aln_file");
	my $line        = <$ALNIN>;
	my $total_chars = 0;
	my $max_chars   = 80;
	while ( $line = <$ALNIN> ) {
		last if $line =~ /^>/;
		chomp $line;
		if ( $max_chars > length($line) ) {
			$line =~ s/-/A/g;
			print $TMPREAD "$line\n";
		} else {
			my $begin = substr( $line, 0, $max_chars - $total_chars );
			$begin =~ s/-/A/g;
			my $end = "-" x ( length($line) - length($begin) );
			print $TMPREAD "$begin$end\n";
		}
		$total_chars += $line;
	}
	return $file.".tmpread.fasta";
}

sub make_ncbi_subtree {
	my %args     = @_;
	my $out_file = $args{out_file} || miss("out_file");
	my $tidref   = $args{taxon_ids} || miss("taxon_ids");

	my @taxonids = @$tidref;

	my ( $nimref, $inmref ) = Phylosift::Summarize::read_ncbi_taxon_name_map();
	my %nameidmap = %$nimref;
	my %idnamemap = %$inmref;
	my $parent    = Phylosift::Summarize::read_ncbi_taxonomy_structure();
	debug "ncbi tree has lots of nodes\n";
	my $merged = Phylosift::Summarize::read_merged_nodes();
	debug "Read a bunch of merged nodes\n";

	my %tidnodes;
	my $phylotree  = Bio::Phylo::Forest::Tree->new();
	my $ncbi_count = 0;
	my $count      = 0;

	foreach my $tid (@taxonids) {
		debug "COULD NOT FIND : ".$tid."\n" unless ( defined( $idnamemap{$tid} ) );
		next unless ( defined( $merged->{$tid} ) || defined( $idnamemap{$tid} ) );    # ensure the id actually exists in NCBI's db
		next if ( $tid eq "" );
		$ncbi_count++;
		my @children;
		while ( defined($tid) ) {

			# check if we've already seen this one
			last if ( defined( $tidnodes{$tid} ) );

			#debug "ADDING $tid\n";
			# process any merging that may have been done
			my @mtid;
			while ( defined( $merged->{$tid} ) ) {
				$tid = $merged->{$tid};

				#				push( @mtid, $tid );
			}
			push( @mtid, $tid );

			# create a new node & add to tree
			my $parentid;
			$parentid = $parent->{$tid}->[0] unless $tid == 1;
			if ( !defined($parentid) && $tid != 1 ) {
				print STDERR "Could not find parent for $tid\n";
				exit;
			}
			my $newnode;
			my @new_children;
			foreach my $mnode (@mtid) {
				if ( defined($parentid) && defined( $tidnodes{$parentid} ) ) {
					$newnode = Bio::Phylo::Forest::Node->new( -parent => $tidnodes{$parentid}, -name => $mnode );
				} else {

					$newnode = Bio::Phylo::Forest::Node->new( -name => $mnode );
				}
				$count++;
				$tidnodes{$mnode} = $newnode;

				# add all children to the new node
				foreach my $child (@children) {
					$newnode->set_child($child);
				}
				$phylotree->insert($newnode);
				push( @new_children, $newnode );
			}

			# continue traversal toward root
			$tid      = $parentid;
			@children = @new_children;
		}
	}

	# if there's something in the tree, write it out
	debug "NCBI COUNT : $ncbi_count\n";
	debug "Making $count Nodes\n";
	if ( $ncbi_count > 0 ) {
		my $TREEOUT = ps_open(">$out_file");
		print $TREEOUT $phylotree->to_newick( "-nodelabels" => 1 );
		close $TREEOUT;
		return 1;    # success
	}
	debug "NOTHING IN THE TREE\n";
	return 0;        # no tree!
}

=head2 generate_id_to_taxonid_map

reads in a hash of unique IDs and a mapping file and writes the <uniqueID><tab><TaxonID> to the output file
the mapping file should be in the following format : <GI><space><TaxonID>
Also returns a list of Taxon IDs
=cut

sub generate_id_to_taxonid_map {
	my %args     = @_;
	my $self     = $args{self} || miss("PS object");
	my $list_ref = $args{id_hash_ref} || miss("Hash for seqID to unique IDs");
	my $map_file = $args{map_file} || miss("seqID to TaxonID file");
	my $out_file = $args{out_file} || miss("ID_taxon mapping filename");
	my $marker   = $args{marker} || miss("Marker name");
	my %map_hash = %{$list_ref};
	my $FHOUT    = Phylosift::Utilities::ps_open( ">".$out_file );
	debug("MAPFILE : $map_file\n");
	my $FHIN = Phylosift::Utilities::ps_open($map_file);

	my @test_array  = keys(%map_hash);
	my %return_hash = ();

	while (<$FHIN>) {
		chomp($_);
		my @line = split( /\t/, $_ );

		#		$_ =~ m/^(\S+)\s+(\S+)$/;
		my $id      = $line[0];
		my $taxonid = $line[1];
		if ( exists $map_hash{$id} ) {
			print $FHOUT $marker."\t".$taxonid."\t".$map_hash{$id}."\n";
			if ( exists $return_hash{$taxonid} ) {
				$return_hash{$taxonid}++;
			} else {
				$return_hash{$taxonid} = 1;
			}
		}

		#no need to print anything otherwise
	}

	#removing any taxon that was seen more than once
	if ($Phylosift::Settings::remove_dup) {
		debug "Removing duplicate Taxons\n";
		foreach my $key ( sort ( keys(%return_hash) ) ) {
			if ( $return_hash{$key} > 1 ) {
				delete( $return_hash{$key} );
			}
		}
	}
	debug "Return_hash has : ".scalar( keys(%return_hash) )." taxons\n";
	return %return_hash;
}

=head2 tree_mangler

mangles a Jplace file according the regex and prints out a tree

=cut

sub tree_mangler {
	my %args     = @_;
	my $fileIN   = $args{in} || miss("Input file");
	my $fileOUT  = $args{out} || miss("Output file");
	my $FH       = Phylosift::Utilities::ps_open($fileIN);
	my $FHOUT    = Phylosift::Utilities::ps_open( ">".$fileOUT );
	my @treedata = <$FH>;
	close($FH);
	my $json_data = decode_json( join( "", @treedata ) );

	# parse the tree
	my $tree_string = $json_data->{tree};
	$tree_string =~ s/:([^:\)\(]+?)\{(\d+?)\}/\{$2\}:$1/g;
	print $FHOUT $tree_string;
}

=head2 generate_hmm

input: alignment_file, target_directory
generates a HMM profile from an alignment in FASTA format (arg) using hmmbuild.  The hmm is placed in the target_directory

=cut

sub generate_hmm {
	my %args      = @_;
	my $file_name = $args{file_name} || miss("file_name");
	my $hmm_name  = $args{hmm_name} || miss("hmm_name");
	`$Phylosift::Settings::hmmbuild --informat afa "$hmm_name" "$file_name"`;
}

=head2 hmmalign_to_model

input : hmm_profile,sequence_file,target_dir
Aligns sequences to an HMM model and outputs an alignment
=cut

sub hmmalign_to_model {
	my %args          = @_;
	my $hmm_profile   = $args{hmm_profile} || miss("hmm_profile");
	my $sequence_file = $args{sequence_file} || miss("sequence_file");
	my $target_dir    = $args{target_dir} || miss("target_dir");
	my $ref_ali       = $args{reference_alignment} || miss("reference_alignment");
	my $seq_count     = $args{sequence_count} || miss("sequence_count");
	my ( $core_name, $path, $ext ) = fileparse( $sequence_file, qr/\.[^.]*$/ );

	my $ALNOUT = ps_open(">$target_dir/$core_name.aln");
	my $ALNIN  = ps_open("$Phylosift::Settings::hmmalign --mapali $ref_ali --trim --outformat afa $hmm_profile $sequence_file |");
	my $s      = 0;
	while ( my $line = <$ALNIN> ) {
		if ( $line =~ /^>/ ) {
			$s++;
			last if $s > $seq_count;
		}
		print $ALNOUT $line;
	}
	close $ALNOUT;
	return "$target_dir/$core_name.aln";
}

=head2 mask_alignment

input : aln_file
Masks the unaligned columns out of an alignment file. Removes ( and ) from the sequence names 
Also removes duplicate IDs
=cut

sub mask_and_clean_alignment {
	my %args        = @_;
	my $aln_file    = $args{alignment_file} || miss("alignment_file");
	my $output_file = $args{output_file} || miss("output_file");
	my %id_map;    # will store a map of unique IDs to sequence names
	debug "Using $aln_file\n";

	#my $in          = Phylosift::Utilities::open_SeqIO_object( file => $aln_file );
	my $IN          = ps_open("$aln_file");       #can't use Bio::SeqIO because the full headers are not being kept as IDs
	my %s           = ();                         #hash remembering the IDs already printed
	my $FILEOUT     = ps_open(">$output_file");
	my $seq_counter = 0;
	my $current_id;
	my $current_seq = "";
	while (<$IN>) {

		#	while ( my $seq_object = $in->next_seq() ) {
		chomp($_);
		if ( $_ =~ m/^>(.*)$/ ) {
			if ( defined $current_id ) {
				my $unique_id = sprintf( "%09d", $seq_counter++ );
				$id_map{$current_id} = $unique_id;

				#$id  =~ s/\(\)//g;                                                       #removes ( and ) from the header lines
				$current_seq =~ s/[a-z]//g;    # lowercase chars didnt align to model
				$current_seq =~ s/\.//g;       # shouldnt be any dots
				print $FILEOUT ">".$unique_id."\n".$current_seq."\n";
				$current_seq = "";             #reset the sequence to empty.
			}
			$current_id = $1;
		} else {
			$current_seq .= $_;
		}
	}

	#    close(FILEIN);
	close($FILEOUT);
	return %id_map;
}

=head2 generate_fasttree

input: alignment_file,target_directory
generates a tree using fasttree and write the output along with the log/info files to the target directory.

=cut

sub generate_fasttree {
	my %args       = @_;
	my $aln_file   = $args{alignment_file} || miss("alignment_file");
	my $target_dir = $args{target_directory} || miss("target_directory");
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my $FILEHANDLE = ps_open($aln_file);
	my %type       = %{ Phylosift::Utilities::get_sequence_input_type($FILEHANDLE) };
	close($FILEHANDLE);
	my $extra_args = $type{seqtype} eq "dna" ? " -nt -gtr " : "";
	debug "DNA alignment detected\n" if ( $type{seqtype} eq "dna" );
	system("$Phylosift::Settings::fasttree $extra_args -log \"$target_dir/$core.log\" \"$aln_file\" > \"$target_dir/$core.tree\" ");
	return ( "$target_dir/$core.tree", "$target_dir/$core.log" );
}

=head2 get_representatives_from_tree

input : tree file in newick format, directory to write the output to, pruning threshold
uses the PDA program to prune a tree to get representative sequences

=cut

sub get_representatives_from_tree {
	my %args       = @_;
	my $tree_file  = $args{tree} || miss("tree");
	my $target_dir = $args{target_directory} || miss("target_directory");
	my $reps_pd    = $args{reps_pd} || miss("reps_pd");
	my ( $core, $path, $ext ) = fileparse( $tree_file, qr/\.[^.]*$/ );

	#get the number of taxa in the tree
	my $taxa_count = 0;
	my $input_tree = new Bio::TreeIO( -file => $tree_file, -format => "newick" );
	while ( my $tree = $input_tree->next_tree ) {
		for my $node ( $tree->get_nodes ) {
			if ( $node->is_Leaf ) {
				$taxa_count++;
			}
		}
	}

	#pda doesn't seem to want to run if $taxa_count is the number of leaves. Decrementing to let pda do the search.
	$taxa_count--;
	my $pda_cmd = "cd \"$target_dir\";$Phylosift::Settings::pda -g -k $taxa_count -minlen $reps_pd \"$tree_file\" \"$target_dir/$core.pda\"";
	`$pda_cmd`;
	return "$target_dir/$core.pda";
}

=head2 get_fasta_from_pda_representatives 

input pda file and reference fasta file
reads the selected representatives from the pda file and prints the sequences to a new fasta file

=cut

sub get_fasta_from_pda_representatives {
	my %args            = @_;
	my $pda_file        = $args{pda_file} || miss("pda_file");
	my $target_dir      = $args{target_dir} || miss("target_dir");
	my $reference_fasta = $args{fasta_reference} || miss("fasta_reference");
	my $id_map_ref      = $args{id_map} || miss("id_map");
	my $core            = $args{core} || miss("core");
	my %id_map          = %{$id_map_ref};

	#reading the pda file to get the representative IDs
	my $REPSIN        = Phylosift::Utilities::ps_open($pda_file);
	my $taxa_number   = 0;
	my %selected_taxa = ();
	while (<$REPSIN>) {
		chomp($_);
		if ( $_ =~ m/optimal PD set has (\d+) taxa/ ) {
			$taxa_number = $1;
		} elsif ( $_ =~ m/Corresponding sub-tree:/ ) {
			last;
		} elsif ( $taxa_number != 0
				  && scalar( keys(%selected_taxa) ) < $taxa_number )
		{
			$_ =~ m/^(\S+)$/;
			$selected_taxa{$1} = 1;
		}
	}
	close($REPSIN);
	my @lect = keys(%selected_taxa);

	#reading the reference sequences and printing the selected representatives using BioPerl
	my $reference_seqs = Phylosift::Utilities::open_SeqIO_object( file   => $reference_fasta,
																  format => "FASTA" );
	my $representatives_fasta = Phylosift::Utilities::open_SeqIO_object( file   => ">$target_dir/$core.rep",
																		 format => "FASTA" );
	while ( my $ref_seq = $reference_seqs->next_seq ) {
		if ( exists $selected_taxa{ $id_map{ $ref_seq->id } } ) {
			$representatives_fasta->write_seq($ref_seq);
		}
	}
	return "$target_dir/$core.rep";
}

sub pd_prune_fasta {
	my %args         = @_;
	my $tre          = $args{tre} || miss("tre");
	my $distance     = $args{distance} || miss("distance");
	my $fasta        = $args{fasta} || miss("fasta");
	my $pruned_fasta = $args{pruned_fasta} || miss("pruned_fasta");

	# no point in pruning something with fewer than three taxa
	my $seq_count = `grep -c ">" "$fasta"`;
	chomp $seq_count;
	if ( $seq_count < 3 ) {
		`cp "$fasta" "$pruned_fasta"`;
		return;
	}

	my $prune_cmd = "$Phylosift::Settings::pda -k 20000 -g -minlen $distance $tre $tre.pruning.log";
	system("$prune_cmd");

	# read the list of taxa to keep
	my $PRUNE  = ps_open("$tre.pruning.log");
	my $intaxa = 0;
	my %keep_taxa;
	while ( my $line = <$PRUNE> ) {
		$intaxa = 1 if ( $line =~ /optimal PD set has/ );
		next unless $intaxa;
		chomp $line;
		$keep_taxa{$line} = 1;
		last if ( length($line) < 2 );    # taxa set ends with empty line
	}

	# create a pruned alignment fasta
	Phylosift::UpdateDB::filter_fasta( input_fasta => $fasta, output_fasta => $pruned_fasta, keep_taxa => \%keep_taxa );
	`rm "$tre.pruning.log"`;
}

=head2 mask_aln

Remove columns from a sequence alignment that contain too many gaps or that are otherwise unreliable

=cut

sub mask_aln {
	my %args       = @_;
	my $infile     = $args{file} || miss("file");
	my $outfile    = $args{output} || miss("output");
	my $gap_cutoff = 10;
	my $cutoff     = 10;
	my $maskcont = martin_mask(
								input_file => $args{file},
								cutoff     => $cutoff,
								opt_g      => $gap_cutoff
	);
	my %maskseq;
	my @ori_order;
	$maskcont =~ s/^>//;
	my @tempmask = split( />/, $maskcont );

	foreach my $tempmask (@tempmask) {
		my @templine = split( /\n/, $tempmask );
		my $this_ID  = shift(@templine);
		my $this_seq = join( '', @templine );
		$this_seq =~ s/\s+//g;
		my @this_seq = split( //, $this_seq );
		if ( $this_ID ne '_mask' ) { push( @ori_order, $this_ID ) }
		$maskseq{$this_ID} = \@this_seq;
	}
	my %trimseq;
	while ( @{ $maskseq{_mask} } ) {
		my $switch = shift( @{ $maskseq{_mask} } );
		foreach my $key ( keys %maskseq ) {
			if ( $key ne '_mask' ) {
				my $aa = shift( @{ $maskseq{$key} } );
				if ( $switch == 1 ) { $trimseq{$key} .= $aa; }
			}
		}
	}

	# writing out a fasta
	my $TRIMOUT = ps_open(">$outfile");
	foreach my $key (@ori_order) {
		print $TRIMOUT ">".$key."\n";
		my $i;
		for ( $i = 0; $i <= length( $trimseq{$key} ); $i += 80 ) {
			my $substr = substr( $trimseq{$key}, $i, 80 );
			print $TRIMOUT $substr."\n";
		}
	}
	close $TRIMOUT;
}

sub martin_mask {
	my %args = @_;
	my ( $input_file, $cutoff, $opt_g ) = ( $args{input_file}, $args{cutoff}, $args{opt_g} );
	if ( !( length($opt_g) > 1 ) ) { $opt_g = 101; }
	my $self = {};
	$self->{inputfile}  = $input_file;
	$self->{cutoff}     = $cutoff;
	$self->{gap_cutoff} = $opt_g;
	$self = &read_matrix( self => $self );
	$self = &read_alignment( self => $self );
	$self = &calculate_score( self => $self );
	$self = &mask( self => $self );
	my $return_mask = &output( self => $self );
	return $return_mask;
}

sub read_matrix {
	my %args = @_;
	my $self = $args{self};
	my @aa;
	my %matrix;
	my ( $matrix_name, $i, $j );
	my %min        = ();
	my $ori_matrix = &ori_matrix();
	my @mx         = split( /\n/, $ori_matrix );
	while (@mx) {
		$_ = shift(@mx);
		if (/amino_acid_order = "(.+)"/) {
			@aa = split //, $1;
		} elsif (/MATRIX (.+)\[/) {
			$matrix_name = $1;
			$j           = 0;
		} elsif (/\d,/) {
			s/(\s|\}|;)//g;
			my @score = split /,/;
			for ( $i = 0; $i <= $#score; $i++ ) {
				$matrix{$matrix_name}{"$aa[$i]$aa[$j]"} = $score[$i] * 100;
				$matrix{$matrix_name}{"$aa[$j]$aa[$i]"} = $score[$i] * 100;
			}
			$j++;
		}
	}
	for my $key ( keys %matrix ) {
		for my $pair ( keys %{ $matrix{$key} } ) {
			if ( !$min{$key} ) {
				$min{$key} = $matrix{$key}{$pair};
			} elsif ( $matrix{$key}{$pair} < $min{$key} ) {
				$min{$key} = $matrix{$key}{$pair};
			}
		}
	}
	foreach my $key ( keys %min ) {
		if ( $min{$key} < 0 ) {
			foreach my $pair ( keys %{ $matrix{$key} } ) {
				$matrix{$key}{$pair} -= $min{$key};
			}
		}
	}
	$self->{matrix} = \%matrix;
	$self->{aa}     = \@aa;
	return $self;
}

sub read_alignment {
	my %args = @_;
	my $self = $args{self};
	my %seq;
	my @order;
	my $id;
	my $input_file = $self->{inputfile};
	my $MASKIN     = ps_open($input_file);
	while (<$MASKIN>) {
		chop;
		if (/%([\S]+)/) {
			$id = $1;
			push( @order, $id );
		} elsif (/>([\S]+)/) {
			$id = $1;
			push( @order, $id );
		} else {
			$_ =~ s/\s+//g;
			$_ =~ s/\./-/g;    #AED: treat . and ? as gaps
			$_ =~ s/\?/-/g;
			$seq{$id} .= uc($_);
		}
	}
	close($MASKIN);
	$self->{seq}   = \%seq;
	$self->{order} = \@order;
	return $self;
}

sub calculate_score {
	my %args = @_;
	my $self = $args{self};
	my %seq  = %{ $self->{seq} };
	my $seqlength;
	my %seqArray;
	my %matrix = %{ $self->{matrix} };
	my %local_score;
	my %column_score;
	my @aa = @{ $self->{aa} };

	for ( keys %seq ) {
		$seqArray{$_} = [ split //, $seq{$_} ];
		$seqlength = length( $seq{$_} );
	}
	if ( $seqlength <= 20 ) {
		carp "Alignment is too short\n";
	}
	my $numseq = scalar keys %seq;
	for my $i ( 0 .. $seqlength - 1 ) {
		my ( %freq, %profile, %dist, %seqVector ) = ();
		my $dist_mean = 0;
		my $number    = 0;
		my $dist_median;
		my $col;
		my $number_gap;
		foreach my $sequence ( keys %seqArray ) {
			$col .= $seqArray{$sequence}[$i];
			if ( $seqArray{$sequence}[$i] ne "-" ) {
				$freq{ $seqArray{$sequence}[$i] }++;
				$number++;
			} else {
				$number_gap++;
			}
		}
		for my $aa1 (@aa) {
			for my $aps (@aa) {
				$profile{$aa1} += $freq{$aps} * $matrix{"gon250mt"}{"$aa1$aps"}
				  if exists $freq{$aps};
			}

			#			$profile{$aa1} /= $number;
			$profile{$aa1} /= $numseq;
		}
		for my $sequence ( keys %seqArray ) {
			my $c = $seqArray{$sequence}[$i];
			if ( $c ne '-' ) {
				for (@aa) {
					$seqVector{$_} = $matrix{"gon250mt"}{"$c$_"};
				}
				for (@aa) {
					my $diff = $profile{$_} - $seqVector{$_};
					$diff /= 1000;
					$dist{$sequence} += $diff * $diff;
				}
				$dist{$sequence} = sqrt( $dist{$sequence} );
				$dist_mean += $dist{$sequence};

				#			$number ++;
			}
		}
		if ($number) {
			$dist_mean /= $number;
		} else {
			$dist_mean = 0;
		}
		my @sort = sort { $a <=> $b } ( values %dist );
		my $t = $number % 2;
		if ( $t == 0 ) {
			$dist_median = ( $sort[ $number / 2 - 1 ] + $sort[ $number / 2 ] ) / 2;
		} else {
			$dist_median = $sort[ ( $number - 1 ) / 2 ];
		}

		#		$column_score2{$i} = exp (-$dist_mean/4) * 100 *$number/$numseq;
		#		$column_score{$i} = exp (-$dist_median/1.82) * 100 * ($number/$numseq);   # 1.82 make total random alignment score 1.0
		$column_score{$i} = exp( -$dist_median / 3 ) * 100 * ( $number / $numseq );
		my $gap_percent = $number_gap / $numseq * 100;
		my $gap_cutoff  = $self->{gap_cutoff};
		if ( ( $gap_cutoff > 0.00001 ) && ( $gap_percent > $gap_cutoff ) ) {
			$self->{rid_gap}->{$i} = 1;
		}
	}
	$self->{column_score} = \%column_score;
	$self->{local_score}  = \%local_score;
##$self->{column_score2}=\%column_score2;
	$self->{seqlength} = $seqlength;
	$self->{seqArray}  = \%seqArray;
	return $self;
}

sub mask {
	my %args         = @_;
	my $self         = $args{self};
	my %matrix       = %{ $self->{matrix} };
	my %seq          = %{ $self->{seq} };
	my %seqArray     = %{ $self->{seqArray} };
	my %column_score = %{ $self->{column_score} };
	my %local_score  = %{ $self->{local_score} };
	my $mask;
	my $seqlength = $self->{seqlength};

	for my $i ( 0 .. $seqlength - 1 ) {
		if ( $i <= 2 or $i >= $seqlength - 3 ) {
			$local_score{$i} = $column_score{$i};
		} else {
			$local_score{$i} =
			  ( $column_score{ $i - 3 } +
				2 * $column_score{ $i - 2 } +
				3 * $column_score{ $i - 1 } +
				4 * $column_score{$i} +
				3 * $column_score{ $i + 1 } +
				2 * $column_score{ $i + 2 } +
				1 * $column_score{ $i + 3 } ) / 16;
		}
		if ( $column_score{$i} == 0 ) { $local_score{$i} = 0; }
		elsif ( $local_score{$i} / $column_score{$i} > 3 ) {
			my $score_l = $column_score{ $i - 3 } + $column_score{ $i - 2 } + $column_score{ $i - 1 };
			my $score_r = $column_score{ $i + 1 } + $column_score{ $i + 2 } + $column_score{ $i + 3 };
			if (    ($score_r)
				 && ( $score_l / $score_r > 20 or $score_l / $score_r < 0.05 ) )
			{
				$local_score{$i} = $column_score{$i};
			}
		}
###########################################cutoff #######################################
		if (    ( $local_score{$i} >= $cutoff )
			 && ( !( $self->{rid_gap}->{$i} ) ) )
		{
			$mask .= "1";
		} else {
			$mask .= "0";
		}

		#		print $i+1,"\t";
		#		printf "%0.1f\t%0.1f\n", $column_score{$i},$local_score{$i};
	}
	$self->{mask} = $mask;
	return $self;
}

sub output {
	my %args        = @_;
	my $self        = $args{self};
	my %seq         = %{ $self->{seq} };
	my $mask        = $self->{mask};
	my @order       = @{ $self->{order} };
	my $return_mask = '';
	$seq{"_mask"} = $mask;
	push( @order, '_mask' );

	foreach my $key (@order) {
		$return_mask .= ">$key\n";
		for ( my $i = 0; $i < length( $seq{$key} ); $i += 60 ) {
			$return_mask .= substr( $seq{$key}, $i, 60 )."\n";
		}
	}
	return $return_mask;
}

sub ori_matrix {
	my $matrix = "amino_acid_order = \"ABCDEFGHIKLMNPQRSTVWXYZ\";

MATRIX gon250mt[]={
66,
35,50,
38,35,127,
33,50,13,82,
35,50,15,68,74,
19,35,29,5,9,97,
54,35,21,35,29,0,94,
29,35,26,37,37,34,25,90,
29,35,27,9,17,41,5,20,77,
32,35,16,38,43,13,27,54,21,72,
27,35,25,8,16,48,5,22,69,21,77,
30,35,29,15,21,45,11,26,67,25,69,79,
33,50,23,65,56,14,37,43,16,40,15,20,76,
37,35,14,30,31,9,24,27,17,31,19,19,29,101,
33,50,19,56,62,17,28,43,22,45,24,28,55,33,68,
31,35,20,33,37,13,28,54,19,68,20,23,37,29,45,82,
58,35,35,38,36,16,53,33,23,35,21,25,41,37,36,33,65,
54,35,31,35,34,20,43,33,31,35,26,31,38,35,35,33,60,67,
35,35,35,15,22,35,13,21,71,23,62,61,20,23,25,21,28,35,73,
11,35,28,0,6,74,8,29,23,11,30,28,11,1,17,24,13,11,17,145,
35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,50,
20,35,31,16,17,84,8,49,30,21,35,33,25,14,23,23,22,22,27,78,35,102,
35,50,35,50,50,35,35,35,35,35,35,35,50,35,50,35,35,35,35,35,35,35,50,";
	return $matrix;
}

=head1 AUTHOR

Martin Wu, C<< <> >>
Dongying Wu, C<< <> >>
Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::MarkerBuild


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Phylosift-Phylosift>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Phylosift-Phylosift>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Phylosift-Phylosift>

=item * Search CPAN

L<http://search.cpan.org/dist/Phylosift-Phylosift/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation.

See http://dev.perl.org/licenses/ for more information.


=cut

1;    # End of Phylosift::MarkerBuild.pm
