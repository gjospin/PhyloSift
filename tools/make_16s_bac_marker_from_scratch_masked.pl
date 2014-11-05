#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use JSON;
use Bio::Phylo::Forest::Tree;
use Bio::Phylo::IO;
use Bio::AlignIO;
#################
#
# make_16s_bac_marker_from_scratch.pl
#
#
# Author : Guillaume Jospin
# Methodology consultant : Aaron Darling
#
# Creates a PhyloSift 16s marker from scratch using bits of codes from the PhyloSift code base
# This is based on the greengenes 99otu tree to find representative sequences
#
# Assume usage of ssu-align http://selab.janelia.org/software/ssu-align/
#            FastTree https://github.com/gjospin/PhyloSift/blob/master/bin/FastTree
#            pda https://github.com/gjospin/PhyloSift/blob/master/bin/pda
#            Readonciler  https://github.com/gjospin/PhyloSift/blob/master/bin/readconciler
#            taxit https://github.com/fhcrc/taxtastic 
#            stockholm2fasta.pl    https://raw.githubusercontent.com/ihh/dart/master/perl/stockholm2fasta.pl
#            hmmalign
#
#
#################

my $debuglevel = 0; # 


# hard coded paths = BAD.
# Need to modify this later for dynamic use
my $taxonomy_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/gg_13_5_otus/taxonomy/99_otu_taxonomy.txt";
my $gg_tree_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/gg_13_5_otus/trees/gg_13_5_otus_99_annotated.tree.gz";
my $gg_accession_mapping_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/gg_13_5_accessions.txt.gz";
my $gg_to_ncbi_mapping_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/taxonIDs/partial_gg_to_ncbi.taxonIDs";
my $raw_sequence_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/gg_13_5.fasta";
my ($marker_name,$domain,$wdir, $leaf_only);
my $ssu_align_model_dir = "/home/gjospin/software/ssu-align-0.1/models/";
# bacteria-0p1.cm OR archaea-0p1.cm

GetOptions ("debug" => \$debuglevel, #flag
	    "domain=s" => \$domain , #string Archaea or Bacteria
	    "marker-name=s" => \$marker_name, #string for the marker name
	    "working-dir=s" => \$wdir, #string where we will write all files and directories to
	    "leaf-only" => \$leaf_only, # how to work the readconciler part
    );

die "Working directory not initialized. Exiting ...\n" unless $wdir;

debug("Creating the working directory\n");
`mkdir -p $wdir` unless -e "$wdir";

# Find list of the taxa from our particular domain
my $count_all = 0;
my $count_full = 0;
debug("Reading the taxonomy $taxonomy_file\n");
my $OUTTAXA = ps_open(">$wdir/gg_all_taxa_$domain.txt");
my $OUTTAXAFULL = ps_open(">$wdir/gg_full_taxa_$domain.txt");
my $INTAXA = ps_open($taxonomy_file);
my $count=0;
while(<$INTAXA>){
    chomp($_);
    if($_ =~ m/^(\d+)\s+k__$domain.*s__\S+$/){
	print $OUTTAXAFULL "$1\n";
	print $OUTTAXA "$1\n";
	$count_all++;
	$count_full++;
    }elsif($_ =~ m/^(\d+)\s+k__$domain/){
	print $OUTTAXA "$1\n";
	$count_all++;
    }
}
close($INTAXA);
close($OUTTAXA);
debug("Done reading $taxonomy_file\n");


# Need to reformat the tree to work with PDA
my $TREEIN = ps_open("zcat $gg_tree_file |");
my $TEMPTREE = ps_open(">$wdir/gg_fixed_tree.tree");
# 1. delete first 4 lines
<$TREEIN>;
<$TREEIN>;
<$TREEIN>;
<$TREEIN>;
# 2. get rid of following chars:
# perl -p -i -e "s/\n//g" gg.broke.tre
# perl -p -i -e "s/;/_/g" gg.broke.tre
# perl -p -i -e "s/ /_/g" gg.broke.tre
# perl -p -i -e "s/'//g" gg.broke.tre 
while(<$TREEIN>){
    $_ =~ s/\n//g;
    $_ =~ s/;/_/g;
    $_ =~ s/ /_/g;
    $_ =~ s/'//g;
    $_ =~ s/\)[a-z]__[^:]{20,}:/\):/g;
    print $TEMPTREE $_;
}
print $TEMPTREE ";";
close($TREEIN);
close($TEMPTREE);
# Cut the tree to only include the taxa from our domain.

unless(-e "$wdir/gg_all_$domain.pda"){
    # Arbitrary number used 0.7 of the total number
    # nothing worked between 140k and 196k.
    my $temp_count = $count_all * 0.7;
    debug("Running pda -i $wdir/gg_all_taxa_$domain.txt -k $temp_count $wdir/gg_fixed_tree.tree $wdir/gg_all_$domain.pda\n");
    `pda -i $wdir/gg_all_taxa_$domain.txt -k $temp_count $wdir/gg_fixed_tree.tree $wdir/gg_all_$domain.pda`;
}

# Extract the tree from the pda output file to
my %selected_taxa=();
unless(-e "$wdir/gg_all_$domain.tree"){
    my $PDAIN        = ps_open("$wdir/gg_all_$domain.pda");
    my $TREEOUT = ps_open(">$wdir/gg_all_$domain.tree");
    my $taxa_number   = 0;
    my $flag = 0;
    while (<$PDAIN>) {
	chomp($_);
	if ( $_ =~ m/Corresponding sub-tree:/ ) {
	    my $tree_string = <$PDAIN>;
	    print $TREEOUT $tree_string;
	}
    }
    close($PDAIN);
    close($TREEOUT);
}

# Extract the taxa from the new tree that are in the list of full taxa and whatever else to go up to our magic number of 25k
unless(-e "$wdir/gg_fulltaxa_$domain.pda"){
    my $temp_count = 25000;
    debug("Running pda -i $wdir/gg_full_taxa_$domain.txt -k $temp_count $wdir/gg_all_$domain.tree $wdir/gg_fulltaxa_$domain.pda\n");
    `pda -i $wdir/gg_full_taxa_$domain.txt -k $temp_count $wdir/gg_all_$domain.tree $wdir/gg_fulltaxa_$domain.pda`;
}

# create an unaligned fasta file from the raw sequences that match the pda picks.
unless(-e "$wdir/gg_25k_unaligned_$domain.fasta" && 0){
    debug("Generating raw files from pda representatives\n");
    my $PDAIN        = ps_open("$wdir/gg_fulltaxa_$domain.pda");
    my $SEQOUT = open_sequence_file(file=>">$wdir/gg_25k_unaligned_$domain.fasta");
    my $SEQIN = open_sequence_file(file=>$raw_sequence_file);
    my $taxa_number   = 0;
    my $flag = 0;
    while (<$PDAIN>) {
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
    close($PDAIN);
    my $print_flag=0;
    while(<$SEQIN>){
	if($_ =~ m/^>(\d+)\n$/){
	    $print_flag=1 if exists $selected_taxa{$1};
	    $print_flag=0 unless exists $selected_taxa{$1};
	}
	print $SEQOUT $_ if $print_flag;
    }
    close($SEQIN);
    close($SEQOUT);
}
debug("Found ".scalar(keys(%selected_taxa))." from the PDA picks\n");
# Align the raw sequences using ssu-align.
unless(-e "$wdir/gg_25k_ssualign_$domain/gg_25k_ssualign_$domain.bacteria.stk"){
#    debug("Running ssu-align -m $ssu_align_model_dir/bacteria-0p1.cm $wdir/gg_25k_unaligned_$domain.fasta $wdir/gg_25k_ssualign_$domain\n");
#    `ssu-align -m $ssu_align_model_dir/bacteria-0p1.cm $wdir/gg_25k_unaligned_$domain.fasta $wdir/gg_25k_ssualign_$domain`;
}


# transform the stockholm alignment into a fasta alignment so FastTree can read it.
unless(-e "$wdir/gg_25k_$domain.aln"){
    debug("Running stockholm2fasta.pl -g -c 80 < $wdir/gg_25k_ssualign_$domain/gg_25k_ssualign_$domain.bacteria.stk > $wdir/gg_25k_$domain.aln\n");
    `perl stockholm2fasta.pl -g -c 80 < $wdir/gg_25k_ssualign_$domain/gg_25k_ssualign_$domain.bacteria.stk > $wdir/gg_25k_$domain.aln`;
}

#changing . into - for pplacer usage
#perl -pi -e 's/\./-/g' gg_13_5_25k_ssualign/gg_13_5_25k_ssualign.bacteria.aln 
#perl -pi -e 's/U/T/g' 16s_test_ref_package/gg_13_5_25k_ssualign.bacteria.aln
#perl -pi -e 's/u/t/g' 16s_test_ref_package/gg_13_5_25k_ssualign.bacteria.aln
unless(-e "$wdir/gg_25k_$domain.clean"){
#    my $seqOUT = Bio::AlignIO->new(-file=>">$wdir/gg_25k_$domain.clean", -format=>"fasta");
#    my $seqIN = Bio::AlignIO->new(-file=>"$wdir/gg_25k_$domain.aln", -format=>"fasta");
#    while ( my $aln = $seqIN->next_aln() ) {
#	while(my $seq = $aln->next_seq()){
#	    $seq->seq =~ s/U/T/g;
#	    $seq->seq =~ s/[a-z]//g;
#	    $seq->seq =~ s/\.//g;
#	}
#    }
#    while ( my $aln = $seqIN->next_aln() ) {    
#	$seqOUT->write_aln($aln);
#    }
    my $ALNIN = open_sequence_file(file=>"$wdir/gg_25k_$domain.aln");
    my $ALNOUT = open_sequence_file(file=>">$wdir/gg_25k_$domain.clean");
    my $curr_seq = "";
    my $curr_id = ""; 
    while(<$ALNIN>){
	chomp($_);
	if($_=~ m/^>(\S+)/){
	    if($curr_seq ne ""){
		$curr_seq =~ s/\.//g;
		$curr_seq =~ s/U/T/g;
		$curr_seq =~ s/[a-z]//g;
		print $ALNOUT ">$curr_id\n$curr_seq\n";
		$curr_seq ="";
	    }
	    $curr_id = $1;
	}else{
	    $curr_seq .= $_;
	}

    }
    # Don't forget the last sequence
    $curr_seq =~ s/\.//g;
    $curr_seq =~ s/U/T/g;
    $curr_seq =~ s/[a-z]//g;
    print $ALNOUT ">$curr_id\n$curr_seq\n";

    close($ALNIN);
    close($ALNOUT);
}

# Create a tree using FastTree
unless(-e "$wdir/gg_25k_$domain.tree" && -e "$wdir/gg_25k_$domain.log"){
    debug("Running FastTree -nt -gtr -log $wdir/gg_25k_$domain.log < $wdir/gg_25k_$domain.clean > $wdir/gg_25k_$domain.tree\n");
    `FastTree -nt -gtr -log $wdir/gg_25k_$domain.log < $wdir/gg_25k_$domain.clean > $wdir/gg_25k_$domain.tree`;
}

unless(-e "$wdir/test_placement.fasta"){
    # need to extract a test placement read
    my $grep = `grep -n '>' $wdir/gg_25k_$domain.clean | head -n 2`;
    $grep=~ m/>\d+\n(\d+):/;
    my $next_header_line = $1;
    debug("Next header line = $next_header_line\n");
# CHANGE HEADER SO THE ID is differente from REFERENCES
    my $SEQIN = open_sequence_file(file=>"$wdir/gg_25k_$domain.clean");
    my $SEQOUT =open_sequence_file(file=>">$wdir/test_placement.fasta");
    my $line_count=2;
    <$SEQIN>;
    print $SEQOUT ">test_placement\n";
    while(<$SEQIN>){
	print $SEQOUT $_;
	last if $line_count > ($next_header_line-2);
        $line_count++;
    }
    close($SEQIN);
    close($SEQOUT);
}

# Create a package for pplacer (also renaming files for PhyloSift marker compliance).
`cp $wdir/gg_25k_$domain.log $wdir/$marker_name.log`;
`cp $wdir/gg_25k_$domain.tree $wdir/$marker_name.tree`;
`cp $wdir/gg_25k_$domain.clean $wdir/$marker_name.clean`;

debug("Running : taxit create -P $wdir/$marker_name -s $wdir/gg_25k_$domain.log -t $wdir/gg_25k_$domain.tree -f $wdir/gg_25k_$domain.clean -l \"$marker_name\"\n");
my $taxit_cmd = "taxit create -P $wdir/$marker_name -s $wdir/$marker_name.log -t $wdir/$marker_name.tree -f $wdir/$marker_name.clean -l \"$marker_name\"";
`$taxit_cmd`;
#place using the newly create pplacer package
unless(-e "$wdir/test_placement.jplace"){
    debug("Running cd $wdir; pplacer -c $marker_name test_placement.fasta\n");
    `cd $wdir; pplacer -c $marker_name test_placement.fasta`;
}

# Finding the mapping between gg ID and NCBI taxon ID
# $selected_taxa

# This part is using some precomputed taxonID fetching from NCBI.
# on edhar @ /home/gjospin/misc/Phylosift_rRNA_markers/20140219/taxonIDs/partial_gg_to_ncbi.taxonIDs
#my $gg_accession_mapping_file ="/home/gjospin/misc/Phylosift_rRNA_markers/20140219/gg_13_5_accessions.txt.gz";
#my $gg_to_ncbi_mapping_file = "/home/gjospin/misc/Phylosift_rRNA_markers/20140219/taxonIDs/partial_gg_to_ncbi.taxonIDs";

#  GI to NCBI taxon mapping
my %ncbi_map=();
my $GGACC = ps_open($gg_to_ncbi_mapping_file);
while(<$GGACC>){
    chomp($_);
    my @line = split(/\s+/,$_);
    $ncbi_map{$line[1]} = $line[2];
}
close($GGACC);

# GG_ID to GI
my %gg_lookup = ();
debug("zcat $gg_accession_mapping_file |\n");
#exit;
my $INGG = ps_open("zcat $gg_accession_mapping_file |");
<$INGG>;
while(<$INGG>){
    chomp($_);
    my @line = split(/\t/,$_);
    $gg_lookup{$line[0]} = $line[2];
}
my %leftovers = ();
#my $LEFTOVER = ps_open(">$wdir/gg_25k_leftover.txt");
my $TAXONMAPOUT = ps_open(">$wdir/gg_25k_$domain.taxmap");
debug("Found ".scalar(keys(%selected_taxa))." selected taxa without ncbi taxonIDs\n");
foreach my $gg_id( keys %selected_taxa){
    if(exists $gg_lookup{$gg_id}){
	if(exists $ncbi_map{$gg_lookup{$gg_id}} && defined $ncbi_map{$gg_lookup{$gg_id}}){
	    print $TAXONMAPOUT "$gg_id\t$ncbi_map{$gg_lookup{$gg_id}}\n";
	    #debug("gg_id $gg_id\ngg_lookup $gg_lookup{$gg_id} \nncbi_map $ncbi_map{$gg_lookup{$gg_id}}\n");
	}else{
	    $leftovers{$gg_id}=$gg_lookup{$gg_id};
	}
    }
}
debug("Found ".scalar(keys(%leftovers))." leftover taxa without ncbi taxonIDs\n");
#exit;
my $div =0;
open(OUT,">$wdir/lookup_acc_left_overs.$div.efetch");
print OUT "for ACC in ";
my $count_map=0;
foreach my $gi(keys %leftovers){
    print OUT $leftovers{$gi}." ";
    $count_map++;
    if($count_map % 5000 == 0){
	print OUT "\ndo\n";
	print OUT "echo -n -e \"\$ACC\\t\"\n";
	print OUT "curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\${ACC}&rettype=fasta&retmode=xml\" |\\\n";
	print OUT "grep TSeq_taxid |\\\n";
	print OUT "cut -d '>' -f 2 |\\\n";
	print OUT "cut -d '<' -f 1 |\\\n";
	print OUT "tr -d \"\\n\"\n";
	print OUT "echo\n";
	print OUT"done\n";
	close (OUT);
	#`chmod 755 $wdir/lookup_acc_left_overs.$div.efetch ; $wdir/lookup_acc_left_overs.$div.efetch > $wdir/lookup_acc_left_overs.$div.taxonIDs`;
	$div++;
	open(OUT , ">$wdir/lookup_acc_left_overs.$div.efetch");
	print OUT "for ACC in ";
    }
}
print OUT "\ndo\n";
print OUT "echo -n -e \"\$ACC\\t\"\n";
print OUT "curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\${ACC}&rettype=fasta&retmode=xml\" |\\\n";
print OUT "grep TSeq_taxid |\\\n";
print OUT "cut -d '>' -f 2 |\\\n";
print OUT "cut -d '<' -f 1 |\\\n";
print OUT "tr -d \"\\n\"\n";
print OUT "echo\n";
print OUT"done\n";
close(OUT);

`chmod 755 $wdir/lookup_acc_left_overs.$div.efetch ; $wdir/lookup_acc_left_overs.$div.efetch > $wdir/lookup_acc_left_overs.$div.taxonIDs`;

#exit;
my %new_ncbi_lookup=();
my @files = <$wdir/lookup_acc_left_overs.*.taxonIDs>;
foreach my $file(@files){
    open(IN,$file);
    while(<IN>){
	chomp($_);
	my @line = split(/\s+/,$_);
	$new_ncbi_lookup{$line[1]} = $line[2];
    }
    close(IN);
}

foreach my $lo (keys %leftovers){
    if(exists $new_ncbi_lookup{$leftovers{$lo}} && defined $new_ncbi_lookup{$leftovers{$lo}}){
	print $TAXONMAPOUT "$lo\t$new_ncbi_lookup{$leftovers{$lo}}\n";
    }
}
close($TAXONMAPOUT);





# Need to prep for readconciler
my $tree_file = "$wdir/test_placement.jplace";
my $mangled = "$wdir/test_placement.mangled";

open(IN,$tree_file);
my @treedata=  <IN>;
close(IN);

open(OUT,">$mangled");
my $json_data = decode_json( join( "", @treedata ) );

# parse the tree
my $tree_string = $json_data->{tree};
$tree_string =~ s/:([^:\)\(]+?)\{(\d+?)\}/\{$2\}:$1/g;
print OUT $tree_string;
close(OUT);


#create NCBI subtree

#get taxon IDs
my @taxonids=();
my $taxon_id_file = "$wdir/gg_25k_$domain.taxmap";
open(INTAXA,$taxon_id_file);
while(<INTAXA>){
    chomp($_);
    my @line = split(/\t/,$_);
    push(@taxonids,$line[1]);
}
my $out_file = "$wdir/ncbi_subtree.tree";
my ( $nimref, $inmref ) = read_ncbi_taxon_name_map();
my %nameidmap = %$nimref;
my %idnamemap = %$inmref;
my $parent    = read_ncbi_taxonomy_structure();
debug("ncbi tree has lots of nodes\n");
my $merged = read_merged_nodes();
debug("Read a bunch of merged nodes\n");

my %tidnodes;
my $phylotree  = Bio::Phylo::Forest::Tree->new();
my $ncbi_count = 0;
my $count_node      = 0;

foreach my $tid (@taxonids) {
    debug("COULD NOT FIND : ".$tid."\n") unless ( defined( $idnamemap{$tid} ) );
    next unless ( defined( $merged->{$tid} ) || defined( $idnamemap{$tid} ) );    # ensure the id actually exists in NCBI's db
    next if ( $tid eq "" );
    $ncbi_count++;
    my @children;
    while ( defined($tid) ) {

	# check if we've already seen this one
	last if ( defined( $tidnodes{$tid} ) );

	# process any merging that may have been done
	my @mtid;
	while ( defined( $merged->{$tid} ) ) {
	    $tid = $merged->{$tid};

	        #push( @mtid, $tid );
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
	    $count_node++;
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
debug("NCBI COUNT : $ncbi_count\n");
debug("Making $count_node Nodes\n");
if ( $ncbi_count > 0 ) {
    open(my $TREEOUT,">$out_file");
    print $TREEOUT $phylotree->to_newick( "-nodelabels" => 1 );
    close $TREEOUT;
#    return 1;    # success
}
debug("Done preparing files for Readconciler\n");

# Run Readconciler  
# Too many ambiguous reads  - empty output
#unless(-e "$wdir/gg_25k_$domain.taxonmap"){
    #debug("Starting Readconciler\n");
    #`readconciler ncbi_subtree.tree $wdir/test_placement.mangled $wdir/gg_25k_$domain.taxmap $wdir/gg_25k_readconciler_$domain.taxonmap`;
    #debug("Done with Readconciler\n");	
#}


# Leaf only Readconciler
debug("Looking for $wdir/$marker_name/$marker_name.taxonmap\n");
#unless(-e "$wdir/$marker_name/$marker_name.taxonmap" && 0){#|| !$leaf_only){
if(0){
    debug("Starting leaf only readconciling\n");
    my $INJASON = ps_open($tree_file);
    my $INMAP = ps_open($taxon_id_file);
    my $OUTTAXONMAP = ps_open(">$wdir/$marker_name/$marker_name.taxonmap");
    
    my @treedata     = <$INJASON>;
    close($INJASON);
    my %id_map=();
    my %unique_to_original_id=();
    my %map=();

    while(<$INMAP>){
	chomp($_);
	$_ =~ m/^(\S+)\s+(\d+)$/;
	$map{$1} = $2;
    }
    close($INMAP);

    my $json_data             = decode_json( join( "", @treedata ) );
    my $tree_string           = $json_data->{tree};

    $tree_string =~ s/^\s+\"//g;
    $tree_string =~ s/:(.+?)(\{\d+?\})/$2:$1/g;
    $tree_string =~ s/\"\,$//g;
    my $tree = Bio::Phylo::IO->parse( '-string' => $tree_string, '-format' => 'newick' )->first;
    foreach my $node ( @{ $tree->get_entities } ) {
	my $name = $node->get_name;
	next unless $name =~ m/^(\d+)\{(\d+)\}$/;
	#print "|$name|\n";
	if ( exists $map{$1} ) {
	    print $OUTTAXONMAP "$2\t$map{$1}\n";
	}
    }
    close($OUTTAXONMAP);

    
}


#
# Last thing, remove all the files we don't need and prepare PhyloSift marker
# 

# Need to change the ssu alignment to have Ts instead of Us.
#my $seqIN = open_SeqIO_object(file=>"$wdir/gg_25k_ssualign_$domain/gg_25k_ssualign_$domain.stk", format=>"stockholm");
#my $seqOUT = open_SeqIO_object(file=>">$wdir/$marker_name/$marker_name.stk", format=>"stockholm");
#while ( my $aln = $seqIN->next_aln() ) {    
#    $seqOUT->write_aln($aln);
#}

`cp $ssu_align_model_dir/bacteria-0p1.cm $wdir/$marker_name/$marker_name.cm`;
`cp $wdir/gg_25k_ssualign_Bacteria/gg_25k_ssualign_Bacteria.bacteria.stk $wdir/$marker_name/$marker_name.stk`;
debug("cp $wdir/$marker_name/$marker_name.clean $wdir/$marker_name/$marker_name.masked\n");
`cp $wdir/$marker_name.clean $wdir/$marker_name/$marker_name.masked`;
`cp $wdir/gg_25k_unaligned_Bacteria.fasta $wdir/$marker_name/$marker_name.rep`;

# need to make HMM for short read placements
`mkdir $wdir/$marker_name.short_tmp` unless -e "$wdir/$marker_name.short_tmp";
unless(-e "$wdir/$marker_name.short_tmp/$marker_name.hmm"){
    `hmmbuild --dna $wdir/$marker_name.short_tmp/$marker_name.hmm $wdir/$marker_name/$marker_name.stk`;
}

# Trying to build a short marker using hmm align and rebuilding the tree.
unless( -e "$wdir/$marker_name.short_tmp/gg_25k_Bacteria2.hmmalign"){
    `/home/gjospin/software/phylosift_v1.0.1/bin/hmmalign --dna --outformat afa -o $wdir/$marker_name.short_tmp/gg_25k_Bacteria2.hmmalign --mapali $wdir/$marker_name/$marker_name.stk $wdir/$marker_name.short_tmp/$marker_name.hmm $wdir/gg_25k_unaligned_Bacteria.fasta`;
}
# Remove duplicate sequences from the hmmalign output

unless( -e "$wdir/$marker_name.short_tmp/$marker_name.clean"){
    my $ALNIN = open_sequence_file(file=>"$wdir/$marker_name.short_tmp/gg_25k_Bacteria2.hmmalign");
    my $ALNOUT = open_sequence_file(file=>">$wdir/$marker_name.short_tmp/$marker_name.clean");
    
    my %IDs = ();
    my $curr_seq = "";
    my $curr_id = "";
    my $flag = 0; 
    while(<$ALNIN>){
	chomp($_);
	if($_=~ m/^>\d+\/\d-\d+/){
	    $flag =0 ;
	    #skip these
	}elsif($_=~ m/^>(\d+)/){
	    $flag =1 ;
	    if($curr_seq ne ""){
		$curr_seq =~ s/\.//g;
#           $curr_seq =~ s/U/T/g;
		$curr_seq =~ s/[a-z]//g;
#           print STDERR "Length seq : ".length($curr_seq)."\n";
		print $ALNOUT ">$curr_id\n$curr_seq\n" unless exists $IDs{$curr_id};
		$IDs{$curr_id}=1;
		$curr_seq ="";
	    }
	    $curr_id = $1;
	}else{
	    $curr_seq .= $_ if $flag;;
	}
    }
    # Don't forget the last sequence
    $curr_seq =~ s/\.//g;
    $curr_seq =~ s/U/T/g;
    $curr_seq =~ s/[a-z]//g;
    print $ALNOUT ">$curr_id\n$curr_seq\n" unless exists $IDs{$curr_id};
    
    close($ALNIN);
    close($ALNOUT);
}

unless(-e "$wdir/$marker_name.short_tmp/$marker_name.log"){
    `/home/gjospin/software/phylosift_v1.0.1/bin/FastTree -log $wdir/$marker_name.short_tmp/$marker_name.log -nt -gtr < $wdir/$marker_name.short_tmp/$marker_name.clean > $wdir/$marker_name.short_tmp/$marker_name.tree_donotuse`;
}
# Create a package for pplacer (also renaming files for PhyloSift marker compliance).



my $taxit_cmd_short = "taxit create -P \"$wdir/$marker_name.short\" -s \"$wdir/$marker_name.short_tmp/$marker_name.log\" -t \"$wdir/$marker_name.tree\" -s \"$wdir/$marker_name.short_tmp/$marker_name.log\" -f \"$wdir/$marker_name.short_tmp/$marker_name.clean\" -l \"$marker_name\"";
print STDERR "Running $taxit_cmd\n";
`$taxit_cmd_short`;
#`taxit create -P $wdir/$marker_name.short -s $wdir/$marker_name.short_tmp/$marker_name.log -t $wdir/$marker_name.tree -s $wdir/$marker_name.short_tmp/$marker_name.log -f $wdir/$marker_name.short_tmp/$marker_name.clean -l \"$marker_name\"`;

`cp $wdir/$marker_name.short_tmp/$marker_name.hmm $wdir/$marker_name.short`;
`cp $wdir/$marker_name/$marker_name.stk $wdir/$marker_name.short`;
`cp $wdir/$marker_name.short/$marker_name.clean $wdir/$marker_name.short/$marker_name.masked`;
`cp $wdir/$marker_name/$marker_name.taxonmap $wdir/$marker_name.short`;
`cp $wdir/$marker_name.short/$marker_name.hmm $wdir/$marker_name`;






################
#
# SUBROUTINES
#
################

sub leaf_only_readconciler {
    my %args         = @_;
    my $mapping_file = $args{mapping};
    my $jplace       = $args{jplace};
    my $output       = $args{output};
    my $id_map_ref   = $args{id_map_ref};
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



sub debug {
    my $msg = shift;
    my $msglevel = shift || 1;
    print $msg if $debuglevel >= $msglevel;
}

sub open_SeqIO_object {
    my %args   = @_;
    my $format = $args{format} || "FASTA";      #default
    my $file   = $args{file};
    my $io_object;
    if ( exists $args{format} ) {
	$format = $args{format};
    }
    if ( $args{file} =~ /\.gz$/ ) {
	$io_object = Bio::SeqIO->new( -file   => "gzip -cd \"$args{file}\" |",
				      -format => $format );
    } elsif ( $args{file} =~ /\.bz2$/ ) {
	$io_object = Bio::SeqIO->new( -file   => "bzcat \"$args{file}\" |",
				      -format => $format );
    } else {
	$io_object = Bio::SeqIO->new( -file => $args{file}, -format => $format );
    }
    return $io_object;
}

sub open_sequence_file {
    my %args = @_;
    my $file = $args{file};
    my $F1IN;
    if ( $file =~ /\.gz$/ ) {
	$F1IN = ps_open("gzip -cd \"$file\" |");
    } elsif ( $file =~ /\.bz2$/ ) {
	$F1IN = ps_open("bzcat \"$file\" |");
    } elsif ( $file eq "STDIN" ) {
	$F1IN = *STDIN;
    } else {
	$F1IN = ps_open($file);
    }
    return $F1IN;
}


sub ps_open {
    my $open_string = shift;
    my $result = open( my $FH, $open_string );
    unless ($result) {
	my $type = "read from";
	$type = "write to"  if $open_string =~ /^>/;
	$type = "append to" if $open_string =~ /^>>/;
	$type = "append to" if $open_string =~ /^\+>/;
	$type = "pipe"      if $open_string =~ /\|$/;
	$Carp::Verbose = 1;
	die "Unable to $type $open_string";
    }
    return $FH;
}



sub stockholm2fasta {
    my %args     = @_;
    my $columns  = $args{columns} || 80;    # number of columns in fasta
    my $gapped   = $args{gapped} || 1;      # should gapped fasta (aligned) be written?
    my $sorted   = $args{sorted} || 1;      # should sequences be sorted?
    my $STREAMIN = $args{in};
    my @seq;
    my @names;
    my $outbuffer = "";
    my $seq_line  = 0;                      # counter for which line we're on in each block.

    while ( my $line = <$STREAMIN> ) {
	if ( $line !~ /\S/ || $line =~ /^\s*#/ || length($line) < 2 ) {
	    $seq_line = 0;
	    next;
	}
	if ( $line =~ /^\s*\/\// ) {
	    $outbuffer .= printseq(
		columns => $columns,
		sorted  => $sorted,
		seq     => \@seq,
		names   => \@names,
		out     => $args{out}
		);
	} else {
	    chomp $line;
	    my ( $name, $seq ) = split /\s+/, $line;
	    $seq =~ s/[\.\-]//g unless $gapped;
	    $seq[$seq_line] .= $seq;
	    $names[$seq_line] = $name;
	    $seq_line++;
	}
    }
    return $outbuffer;
}


sub printseq {
    my %args  = @_;
    my @seq   = @{ $args{seq} };
    my @names = @{ $args{names} };
    my $out   = "";
    for ( my $j = 0; $j < @names; $j++ ) {
	$out .= ">$names[$j]\n";
	for ( my $i = 0; $i < length( $seq[$j] ); $i += $args{columns} ) {
	    $out .= substr( $seq[$j], $i, $args{columns} )."\n";
	}
    }
    return $out;
}


sub read_merged_nodes {
    my %merged;
    print STDERR "Reading merged ncbi nodes\n";
    open(my $MERGED ,"/home/gjospin/share/phylosift/ncbi/merged.dmp");
    while ( my $line = <$MERGED> ) {
	chomp $line;
	my @vals = split( /\s+\|\s*/, $line );
	$merged{ $vals[0] } = $vals[1];
    }
    print STDERR "Done reading merged\n";
    return \%merged;
}

sub read_ncbi_taxon_name_map {
    my(%nameidmap, %idnamemap );
    my $ncbidir = "/home/gjospin/share/phylosift/ncbi";
    open(my $TAXIDS ,"$ncbidir/names.dmp");
    while ( my $line = <$TAXIDS> ) {
	chomp $line;
	if (    ( $line =~ /scientific name/ )
		 || ( $line =~ /synonym/ )
		 || ( $line =~ /misspelling/ ) )
	{
	    my @vals = split( /\s+\|\s+/, $line );
	    $nameidmap{ homogenize_name_ala_dongying( name => $vals[1] ) } = $vals[0];
	    $idnamemap{ $vals[0] } = homogenize_name_ala_dongying( name => $vals[1] )
		if ( $line =~ /scientific name/ );
	}
    }
    return ( \%nameidmap, \%idnamemap );
}

sub read_ncbi_taxonomy_structure {
    my %parent;

    my $ncbidir      = "/home/gjospin/share/phylosift/ncbi";
    open(my $TAXSTRUCTURE,"$ncbidir/nodes.dmp");
    while ( my $line = <$TAXSTRUCTURE> ) {
	chomp $line;
	my @vals = split( /\s+\|\s+/, $line );
	$parent{ $vals[0] } = [ $vals[1], $vals[2] ];
    }
    return \%parent;
}

sub homogenize_name_ala_dongying {
    my %args = @_;
    my $inName = $args{name} || miss("name");
    return "" unless defined($inName);
    $inName =~ s/^\s+//;
    $inName =~ s/\s+$//;
    $inName =~ s/\s+/ /g;
    $inName =~ s/,//g;
    $inName =~ s/[)(]//g;
    $inName =~ s/-/ /g;
    $inName = uc $inName;
    return $inName;
}
