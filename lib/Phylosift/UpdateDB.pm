package Phylosift::UpdateDB;
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Carp;
use Phylosift::Summarize;
use Bio::Phylo::IO qw(parse unparse);

=head1 NAME

Phylosift::UpdateDB - Functionality to download new genomes and update the marker database

=head1 VERSION

Version 0.01

=cut
our $VERSION = '0.01';

=head1 SYNOPSIS

Acquires sequence data from either a remote or local data repository, then adds that data to the current 
set of markers to (hopefully) improve taxonomic resolution.

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 get_ebi_genomes

=cut

sub get_ebi_genomes {
	my %args      = @_;
	my $directory = $args{directory} // miss("directory");
	chdir($directory);
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/organelle.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/virus.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/phage.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaea.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaealvirus.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/bacteria.details.txt" );

	#	get_ebi_from_list("http://www.ebi.ac.uk/genomes/eukaryota.details.txt");
}

sub get_ebi_from_list() {
	my %args     = @_;
	my $list_url = $args{url} // miss("url");    # URL to EBI's table of genome characteristics
	`wget $list_url -O list.txt`;
	open( DETAILS, "list.txt" );
	my $line = <DETAILS>;
	while ( $line = <DETAILS> ) {

		# each line in the list file records a genome and some metadata
		chomp $line;
		my ( $acc, $version, $date, $taxid, $description ) = split( /\t/, $line );
		my $outfile = $acc;
		$outfile =~ s/\./_/g;
		$outfile .= ".$taxid";
		print STDERR "$outfile.fasta\n";

		# if it exists and is newer, don't re-download
		my $download = 1;
		if ( -e "$outfile.embl" ) {
			my $mtime    = ( stat("$outfile.embl") )[9];
			my @timerval = localtime($mtime);
			my $datestr  = ( 1900 + $timerval[5] );
			$datestr .= 0 if $timerval[4] < 9;
			$datestr .= ( $timerval[4] + 1 );
			$datestr .= 0 if $timerval[3] < 9;
			$datestr .= $timerval[3];
			next unless int($datestr) < int($date);
			`rm $outfile.fasta`;
		}

		# either we don't have this one yet, or our version is out of date
		open( WGETTER, "| wget -i - -O $outfile.embl " );
		print WGETTER "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
		close WGETTER;

		#		eval{
		my $seq_in  = Phylosift::Utilities::open_SeqIO_object( file => "$outfile.embl",   format => "embl" );
		my $seq_out = Phylosift::Utilities::open_SeqIO_object( file => ">$outfile.fasta", format => "fasta" );
		while ( my $inseq = $seq_in->next_seq ) {
			$seq_out->write_seq($inseq);
		}

		#			`rm $outfile.embl`;
		#		} or do {
		#			carp "Error processing $outfile.embl\n";
		#		}
	}
	`rm list.txt`;
}

sub get_taxid_from_gbk() {
	my %args = @_;

	# first get the taxon ID
	my $file = $args{file} // miss("file");
	open( GBK, $file );
	my $taxid;
	while ( my $l2 = <GBK> ) {
		if ( $l2 =~ /\/db_xref\=\"taxon\:(\d+)/ ) {
			$taxid = $1;
			last;
		}
	}
	return $taxid;
}

sub get_ncbi_finished_genomes() {
	my %args      = @_;
	my $directory = $args{directory} // miss("directory");
	`mkdir -p $directory`;

	# First download all finished bacterial genomes
	# then for each genome, concatenate all replicons into a FastA file
	# Name the FastA file with the organism and its taxon ID.
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria";

	#	`$ncbi_wget_cmd`;
	open( FINDER, "find ftp.ncbi.nih.gov/genomes/Bacteria -type d |" );
	while ( my $line = <FINDER> ) {
		chomp $line;
		$line =~ /\/Bacteria\/(.+)/;
		my $orgname = $1;
		next if ( $orgname =~ /\// );    # could be a malformed genbank directory
		next unless length($1) > 1;
		my $seq_out;
		my $fasta_name;
		open( LSSER, "ls $line/*gbk |" );

		while ( my $gbk = <LSSER> ) {
			chomp $gbk;
			if ( !defined($fasta_name) ) {
				my $taxid = get_taxid_from_gbk ( file => $gbk );
				$fasta_name = "$orgname.$taxid.fasta";

				# skip this one if it exists and hasn't been updated
				# POSSIBLE BUG: if only part of this organism's genbank record is updated
				if ( -e $fasta_name && ( stat($fasta_name) )[9] > ( stat($gbk) )[9] ) {
					print STDERR "Already have $fasta_name\n";
					last;
				}
				$seq_out = Phylosift::Utilities::open_SeqIO_object( file => ">$fasta_name", format => "FASTA" );
			}
			my $seq_in = Phylosift::Utilities::open_SeqIO_object( file => "$gbk" );
			while ( my $inseq = $seq_in->next_seq ) {
				$seq_out->write_seq($inseq);
			}
		}
	}
}

sub get_ncbi_draft_genomes() {
	my %args      = @_;
	my $directory = $args{directory} // miss("directory");
	`mkdir -p $directory`;
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	$ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=\"*fna.tgz\",\"*fna.[0-9].tgz\" ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	open( FINDER, "find . -name \"*.gbk\" |" );

	while ( my $line = <FINDER> ) {
		chomp $line;
		my $taxid = get_taxid_from_gbk ( file => $line );

		# then unpack all the files and combine them
		my $fasta_out = basename( $line, ".gbk" ) . ".$taxid.fasta";
		my $fna = $line;
		$fna =~ s/\.gbk/\.scaffold\.fna\.tgz/g;
		$fna =~ s/\.scaffold\.fna\.tgz/\.contig\.fna\.tgz/g unless ( -e $fna );
		$fna =~ s/\.contig\.fna\.tgz/\.contig\.fna\.1\.tgz/g unless ( -e $fna );
		unless ( -e $fna ) {
			warn "Missing FastA data for $line\n";
			next;
		} elsif ( $fna =~ /contig\.fna\.1\.tgz/ ) {

			# if the contigs are broken into many files, unpack them all at once
			$fna =~ s/\.contig\.fna\.1\.tgz/\.contig\.fna\.\*\.tgz/g;
		}

		# skip this one if it exists and hasn't been updated
		if ( -e $fasta_out && ( stat($fasta_out) )[9] > ( stat($line) )[9] ) {
			print STDERR "Already have $fasta_out\n";
			next;
		}

		# unpack the nucleotide tarball and cat the scaffolds/contigs.
		my @tarfiles  = `tar xvzf $fna`;
		my $tarstatus = ( $? >> 8 );
		my $catline   = join( " ", @tarfiles );
		$catline =~ s/\n//g;
		`rm -f $fasta_out`;
		`cat $catline >> $fasta_out` if ( $tarstatus == 0 );
		`rm $catline`;
	}
}

sub find_new_genomes() {
	my %args        = @_;
	my $genome_dir  = $args{genome_directory} // miss("genome_directory");
	my $results_dir = $args{results_directory} // miss("results_directory");
	my $files       = $args{files} // miss("files");
	open( FINDER, "find $genome_dir |" );
	while ( my $genome = <FINDER> ) {
		chomp $genome;
		next unless $genome =~ /\.fasta/;
		my $gbase = basename $genome;
		if ( -e "$results_dir/$gbase/alignDir/concat.trim.fasta" ) {
			my $ctime = ( stat("$results_dir/$gbase/alignDir/concat.trim.fasta") )[9];
			my $mtime = ( stat($genome) )[9];
			push( @{$files}, $genome ) if ( $ctime < $mtime );

			#			push( @{$files}, $genome );
		} else {
			push( @{$files}, $genome );
		}
	}
}

# limit on the number of jobs that can be queued to SGE at once
use constant MAX_SGE_JOBS => 5000;

sub qsub_updates() {
	my %args        = @_;
	my $results_dir = $args{results_directory} // miss("results_directory");
	my $files       = $args{files} // miss("files");
	my @jobids;
	`mkdir -p $results_dir`;
	open( PHYLOSIFTSCRIPT, ">/tmp/pssge.sh" );
	print PHYLOSIFTSCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
#\$ -q eisen.q -q all.q
export PATH=\$PATH:/home/koadman/development/PhyloSift/bin/
export PERL5LIB=\$HOME/lib/perl5:/home/koadman/development/PhyloSift/lib
WORKDIR=/state/partition1/koadman/phylosift/\$JOB_ID
mkdir -p \$WORKDIR
CURDIR=`pwd`
cd \$WORKDIR
phylosift search -isolate -besthit \$1
phylosift align -isolate -besthit \$1
rm -rf PS_temp/*/treeDir
rm -rf PS_temp/*/blastDir
rm -rf PS_temp/*/isolates.fasta
rm -rf PS_temp/*/alignDir/PMP*.stk
rm -rf PS_temp/*/alignDir/PMP*.hmm
rm -rf PS_temp/*/alignDir/PMP*hmm.fasta
rm -rf PS_temp/*/alignDir/PMP*.newCandidate

# ensure orderly copying back to shared storage
LOCKFILE=\$CURDIR/copylock
while [ ! `mktemp -q \$LOCKFILE` ]; do
        sleep 10s
done

cp -R PS_temp/* \$CURDIR/

rm \$LOCKFILE   # remove the lock now that copying is done
rm -rf \$WORKDIR
};
	my $job_count = 0;

	foreach my $file ( @{$files} ) {
		$job_count++;
		chdir($results_dir);
		my $job = `qsub /tmp/pssge.sh $file`;
		$job =~ /Your job (\d+) /;
		push( @jobids, $1 );

		# check whether we've hit the limit for queued jobs, and rest if needed
		if ( $job_count == MAX_SGE_JOBS ) {
			wait_for_jobs(@jobids);
			@jobids    = ();
			$job_count = 0;
		}
	}
	wait_for_jobs( job_ids => @jobids );
}

sub wait_for_jobs {
	my %args = @_;
	my $job_ids = $args{job_ids}  // miss("job_ids");

	# wait for all jobs to complete
	while ( my $jobid = $job_ids ) {
		while (1) {
			my $output = `qstat -j $jobid 2>&1`;
			last if $output =~ /Following jobs do not exist/;
			sleep(20);
		}
	}
}

sub collate_markers() {
	my %args        = @_;
	my $results_dir = $args{results_dir} // miss("results_dir");
	my $marker_dir  = $args{marker_dir} // miss("marker_dir");

	# get list of markers
	chdir($marker_dir);
	my @markerlist = Phylosift::Utilities::gather_markers();
	print STDERR "Markers are " . join( " ", @markerlist ) . "\n";

	# get a list of genomes available in the results directory
	# this hopefully means we touch each inode over NFS only once
	# NFS is slooooow...
	print STDERR "Working with " . scalar(@markerlist) . " markers\n";
	print STDERR "Listing all files in results dir\n";
	my @alldata = `find $results_dir -name "*.trim.fasta"`;
	print STDERR "Found " . scalar(@alldata) . " files\n";
	foreach my $marker (@markerlist) {
		my $cat_ch = ">";    # first time through ensures that existing files get clobbered

		# find all alignments with this marker
		my @catfiles = ();
		foreach my $file (@alldata) {
			next unless $file =~ /(.+\.fasta)\/alignDir\/$marker.trim.fasta/;
			my $genome = $1;
			chomp($file);
			push( @catfiles, $file );

			# cat the files into the alignment in batches
			# this cuts down on process spawning and file I/O
			# can't do more than about 4k at once because argument lists get too long
			if ( @catfiles > 200 ) {
				my $catline = join( " ", @catfiles );
				print STDERR "Found 200 files for marker $marker\n";
				`cat $catline $cat_ch $marker_dir/$marker.updated.fasta`;
				$catline =~ s/\.trim\.fasta/\.trim\.fna\.fasta/g;
				if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
					`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
				}
				$catline =~ s/\.trim\.fna\.fasta/\.unmasked/g;
				`cat $catline $cat_ch $marker_dir/$marker.updated.reps`;
				@catfiles = ();
				$cat_ch   = ">>";    # now add to existing files
			}
		}
		if ( @catfiles > 0 ) {
			my $catline = join( " ", @catfiles );
			print STDERR "Last cat for marker $marker\n";
			`cat $catline $cat_ch $marker_dir/$marker.updated.fasta`;
			$catline =~ s/\.trim\.fasta/\.trim\.fna\.fasta/g;
			if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
				`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
			}
			$catline =~ s/\.trim\.fna\.fasta/\.unmasked/g;
			`cat $catline $cat_ch $marker_dir/$marker.updated.reps`;
		}
		fix_names_in_alignment( alignment => "$marker_dir/$marker.updated.fasta" );
		fix_names_in_alignment( alignment => "$marker_dir/$marker.codon.updated.fasta" );
	}
	`rm -f $marker_dir/concat.updated.fasta`;
	`rm -f $marker_dir/concat.codon.updated.fasta`;
	`cat $results_dir/*/alignDir/concat.fasta >> $marker_dir/concat.updated.fasta`;
	`cat $results_dir/*/alignDir/concat-dna.fasta >> $marker_dir/concat.codon.updated.fasta`;
	fix_names_in_alignment( alignment => "$marker_dir/concat.updated.fasta" );
	fix_names_in_alignment( alignment => "$marker_dir/concat.codon.updated.fasta" );
}

sub clean_representatives {
	my %args    = @_;
	my $file    = $args{infile} // miss("infile");
	my $outfile = $args{outfile} // miss("outfile");
	open( INALN,  $file );
	open( OUTALN, ">$outfile" );
	while ( my $line = <INALN> ) {
		unless ( $line =~ /^>/ ) {
			my $first_lower;
			my $last_upper;

			# strip all chars outside the first and last uppercase character -- these are the boundaries of the homologous region
			for ( my $i = 0 ; $i < 2 ; $i++ ) {
				if ( $line =~ /[A-Z]/ ) {
					$line = substr( $line, $-[0] );
				}
				$line = reverse($line);
			}
			$line =~ s/-//g;     # remove gap chars
			$line =~ s/\.//g;    # remove gap chars
			$line .= "\n";       # trailing newline was removed above
		}
		print OUTALN $line;
	}
}

=head2 assign_seqids
 assign 10 digit unique identifiers to each gene
 todo: will need to move into alphabetical space here
=cut

sub assign_seqids($) {
	my %args       = @_;
	my $marker_dir = $args{marker_directory} // miss("marker_directory");

	# get list of markers
	chdir($marker_dir);
	my @markerlist    = Phylosift::Utilities::gather_markers();
	my $aa_counter    = 0;
	my $codon_counter = 0;
	open( my $AA_IDTABLE,    ">gene_ids.aa.txt" );
	open( my $CODON_IDTABLE, ">gene_ids.codon.txt" );
	push( @markerlist, "concat" );
	foreach my $marker (@markerlist) {
	}

	sub assign_seqids_for_marker() {
		my %args         = @_;
		my $marker       = $args{marker};
		my $alignment    = $args{alignment};
		my $idtable      = $args{id_table};
		my $counter      = $args{counter};
		my $existing_ids = $args{existing_ids};
		open( INALN,  $alignment )           || croak "Unable to read $alignment";
		open( OUTALN, ">$alignment.seqids" ) || croak "Unable to write to $alignment.seqids";
		my %mapped_ids;
		my %id_mapping       = ();
		my %id_mapping_codon = ();
		$aa_counter = assign_seqids_for_marker ( marker => $marker, alignment => "$marker.updated.fasta", id_table => $aa_idtable, counter => $aa_counter,
												 existing_ids => \%id_mapping );
		$codon_counter = assign_seqids_for_marker (
													marker       => $marker,
													alignment    => "$marker.codon.updated.fasta",
													id_table     => $codon_idtable,
													counter      => $codon_counter,
													existing_ids => \%id_mapping_codon
		);
		my $printing = 0;

		while ( my $line = <INALN> ) {
			chomp $line;
			if ( $line =~ /^>(.+)/ ) {
				my $header = $1;
				if ( defined( $mapped_ids{$header} ) ) {

					# this one was already included
					$printing = 0;
					next;
				}
				$printing = 1;
				my $countstring = sprintf( "%010u", $counter );
				if ( defined( $existing_ids->{$header} ) ) {

					# this one uses an ID created previously
					$countstring = $existing_ids->{$header};
				} else {

					# record the new ID in the table
					print $IDTABLE "$marker\t$header\t$countstring\n";
				}
				$mapped_ids{$header} = $countstring;
				$line = ">$countstring";
				$counter++;
			}
			print OUTALN $line . "\n" if $printing;
		}
		close INALN;
		close OUTALN;
		`mv $alignment.seqids $alignment`;
		foreach my $key ( keys %mapped_ids ) {
			$existing_ids->{$key} = $mapped_ids{$key};
		}
		return $counter;
	}

	sub read_gene_ids {
		my %args = @_;
		my $file = $args{file} // miss("file");
		open( IDS, $file ) || croak "Unable to read $file\n";
		my %id_to_taxon;
		my %marker_taxon_to_id;
		while ( my $line = <IDS> ) {
			my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
			$id_to_taxon{$uniqueid} = $taxon;
			$marker_taxon_to_id{$marker}{$taxon} = [] unless defined( $marker_taxon_to_id{$marker}{$taxon} );
			push( @{ $marker_taxon_to_id{$marker}{$taxon} }, $uniqueid );
		}
		return ( \%id_to_taxon, \%marker_taxon_to_id );
	}

	sub update_ncbi_taxonomy {
		my %args = @_;
		my $repository = $args{repository} ||  // miss("repository");
		print "Downloading new NCBI taxonomy...\n";
		my $url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
		my $ff = File::Fetch->new( uri => $url );
		`mkdir -p $repository/ncbi`;
		$ff->fetch( to => "$repository/ncbi" );
		system("cd $repository/ncbi/ ; tar xzf taxdump.tar.gz");
		unlink("$repository/ncbi/taxdump.tar.gz");
		unlink("$repository/ncbi/citations.dmp");
		unlink("$repository/ncbi/delnodes.dmp");
		unlink("$repository/ncbi/division.dmp");
		unlink("$repository/ncbi/gc.prt");
		unlink("$repository/ncbi/gencode.dmp");
		unlink("$repository/ncbi/readme.txt");
		system("cd $repository/; tar czf ncbi.tar.gz ncbi");

		# force use of the new NCBI data
		$Phylosift::Utilities::ncbi_dir = "$repository/ncbi/";
	}

	sub read_merged_nodes {
		open( MERGED, "$Phylosift::Utilities::ncbi_dir/merged.dmp" ) || croak "Unable to read $Phylosift::Utilities::ncbi_dir/merged.dmp";
		my %merged;
		while ( my $line = <MERGED> ) {
			chomp $line;
			my @vals = split( /\s+\|\s*/, $line );
			$merged{ $vals[0] } = $vals[1];
		}
		return %merged;
	}

	sub make_ncbi_tree_from_update {
		my %args        = @_;
		my $self        = $args{self} // miss("self");
		my $results_dir = $args{results_dir} // miss("results_dir");
		my $markerdir   = $args{marker_dir} // miss("marker_dir");
		my ( %nameidmap, %idnamemap ) = Phylosift::Summarize::readNcbiTaxonNameMap();
		my %parent = Phylosift::Summarize::readNcbiTaxonomyStructure();
		print STDERR "ncbi tree has " . scalar( keys(%parent) ) . " nodes\n";
		open( AAIDS, "$markerdir/gene_ids.aa.txt" ) || croak "Unable to read $markerdir/gene_ids.aa.txt";
		open( MARKERTAXONMAP, ">$markerdir/marker_taxon_map.updated.txt" ) || croak "Unable to write $markerdir/marker_taxon_map.updated.txt";
		my @taxonids;
		my %merged = read_merged_nodes();
		print "Read " . scalar( keys(%merged) ) . " merged nodes\n";

		while ( my $line = <AAIDS> ) {
			chomp $line;
			my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
			push( @taxonids, $taxon ) if $taxon =~ /^\d+$/;
			print MARKERTAXONMAP "$uniqueid\t$taxon\n";
		}
		close MARKERTAXONMAP;
		my %tidnodes;
		my $phylotree = Bio::Phylo::Forest::Tree->new();
		foreach my $tid (@taxonids) {
			next if ( $tid eq "" );
			my @children;
			while ( $tid != 1 ) {

				# check if we've already seen this one
				last if ( defined( $tidnodes{$tid} ) );

				# process any merging that may have been done
				my @mtid;
				push( @mtid, $tid );
				while ( defined( $merged{$tid} ) ) {
					$tid = $merged{$tid};
					push( @mtid, $tid );
				}

				# create a new node & add to tree
				my $parentid = $parent{$tid}->[0];
				if ( !defined($parentid) ) {
					print STDERR "Could not find parent for $tid\n";
					exit;
				}
				my $newnode;
				my @new_children;
				foreach my $mnode (@mtid) {
					$newnode = Bio::Phylo::Forest::Node->new( -parent => $tidnodes{$parentid}, -name => $mnode ) if defined( $tidnodes{$parentid} );
					$newnode = Bio::Phylo::Forest::Node->new( -name => $mnode ) if !defined( $tidnodes{$parentid} );
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
		open( TREEOUT, ">ncbi_tree.updated.tre" );
		print TREEOUT $phylotree->to_newick( "-nodelabels" => 1 );
		close TREEOUT;
	}

	sub get_marker_name_base {
		my %args = @_;
		my $name = $args{marker} // miss("marker");
		$name .= ".codon"   if $args{dna} // miss("dna");
		$name .= ".updated" if $args{updated} // miss("updated");

		# rna markers don't get pruned
		$name .= ".pruned" if $args{pruned} && !Phylosift::Utilities::is_protein_marker( marker => $name );
		return $name;
	}

	sub get_fasta_filename {
		my %args = @_;
		my $name = get_marker_name_base(%args);
		$name .= ".fasta";
		return $name;
	}

	sub get_taxonmap_filename {
		my %args = @_;
		my $name = get_marker_name_base(%args);
		$name .= ".taxonmap";
		return $name;
	}

	sub get_reps_filename {
		my %args  = @_;
		my $clean = $args{clean} || 0;
		my $name  = get_marker_name_base(%args);
		$name .= ".reps";
		$name .= ".clean" if $clean;
		return $name;
	}

	sub get_fasttree_log_filename {
		my %args = @_;
		my $name = get_marker_name_base(%args);
		$name .= ".fasttree.log";
		return $name;
	}

	sub get_fasttree_tre_filename {
		my %args = @_;
		my $name = get_marker_name_base(%args);
		$name .= ".tre";
		return $name;
	}

	sub build_marker_trees_fasttree() {
		my %args        = @_;
		my $marker_dir  = $args{directory} // miss("directory");
		my $pruned      = $args{pruned} // miss("pruned");
		my $codon_fasta = get_fasta_filename( marker => '$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_fasta    = get_fasta_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );
		my $codon_tre   = get_fasttree_tre_filename( marker => '$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_tre      = get_fasttree_tre_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );
		my $codon_log   = get_fasttree_log_filename( marker => '$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_log      = get_fasttree_log_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );
		open( TREESCRIPT, ">/tmp/ps_tree.sh" );
		print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -log $aa_log $aa_fasta > $aa_tre
};
		`chmod 755 /tmp/ps_tree.sh`;
		close TREESCRIPT;
		open( TREESCRIPT, ">/tmp/ps_tree_codon.sh" );
		print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -nt -gtr -log $codon_log  $codon_fasta > $codon_tre
};
		`chmod 755 /tmp/ps_tree_codon.sh`;
		close TREESCRIPT;
		open( TREESCRIPT, ">/tmp/ps_tree_rna.sh" );
		print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -nt -gtr -log $aa_log $aa_fasta > $aa_tre
};
		`chmod 755 /tmp/ps_tree_rna.sh`;
		chdir($marker_dir);
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		my @jobids;

		foreach my $marker (@markerlist) {
			my $marker_fasta = get_fasta_filename( marker => $marker, dna => 0, updated => 1, pruned => $pruned );
			print "Looking for $marker_fasta\n";
			next unless ( -e $marker_fasta );
			my $tree_script = "/tmp/ps_tree.sh";
			$tree_script = "/tmp/ps_tree_rna.sh" unless Phylosift::Utilities::is_protein_marker( marker => $marker );

			# run fasttree on them
			my $qsub_cmd = "qsub -q all.q -q eisen.q $tree_script $marker";
			my $job      = `$qsub_cmd`;
			$job =~ /Your job (\d+) /;
			push( @jobids, $1 );
			next unless Phylosift::Utilities::is_protein_marker( marker => $marker );

			# run fasttree on codons
			$qsub_cmd = "qsub -q all.q -q eisen.q /tmp/ps_tree_codon.sh $marker";
			$job      = `$qsub_cmd`;
			$job =~ /Your job (\d+) /;
			push( @jobids, $1 );
		}
		wait_for_jobs( job_ids => @jobids );
		`rm ps_tree.sh.*`;
	}

	sub build_marker_trees_raxml() {
		my %args       = @_;
		my $marker_dir = $args{marker_directory} // miss("marker_directory");
		open( TREESCRIPT, ">/tmp/ps_tree.sh" );
		print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -pe threaded 3
#\$ -S /bin/bash

rm RAxML*\$1*
raxmlHPC -m GTRGAMMA -n \$1.codon.updated -s \$1.codon.updated.phy -T 4
raxmlHPC -m PROTGAMMAWAGF -n \$1.updated -s \$1.updated.phy -T 4
mv RAxML_bestTree.\$1.codon.updated \$1.codon.updated.tre
mv RAxML_bestTree.\$1.updated \$1.updated.tre
mv RAxML_info.\$1.codon.updated \$1.codon.updated.RAxML_info
mv RAxML_info.\$1.updated \$1.updated.RAxML_info
rm RAxML*.\$1*
};
		`chmod 755 /tmp/ps_tree.sh`;
		chdir($marker_dir);
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		my @jobids;

		foreach my $marker (@markerlist) {

			# convert to phylip
			next unless ( -e "$marker.updated.fasta" );
			my $in  = Bio::AlignIO->new( -file => "$marker.updated.fasta", -format => "fasta" );
			my $out = Bio::AlignIO->new( -file => ">$marker.updated.phy",  -format => "phylip" );
			while ( my $inseq = $in->next_aln ) {
				$out->write_aln($inseq);
			}
			if ( -e "$marker.codon.updated.fasta" ) {
				fix_names_in_alignment( alignment => "$marker_dir/$marker.codon.updated.fasta" );
				my $in_dna  = Bio::AlignIO->new( -file => "$marker.codon.updated.fasta", -format => "fasta" );
				my $out_dna = Bio::AlignIO->new( -file => ">$marker.codon.updated.phy",  -format => "phylip" );
				while ( my $inseq = $in_dna->next_aln ) {
					$out_dna->write_aln($inseq);
				}
			}

			# run raxml on them
			my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/ps_tree.sh $marker";
			my $job      = `$qsub_cmd`;
			$job =~ /Your job (\d+) /;
			push( @jobids, $1 );
		}
		wait_for_jobs( job_ids => @jobids );
		`rm ps_tree.sh.*`;
	}

	sub fix_names_in_alignment() {
		my %args      = @_;
		my $alignment = $args{alignment} // miss("alignment");
		my %markertaxa;

		# naive means to remove most trivial duplicates. TODO: allow divergent paralogs to remain.
		my $printing = 1;
		open( INALN,  $alignment );
		open( OUTALN, ">$alignment.fixed" );
		while ( my $line = <INALN> ) {
			if ( $line =~ /^>(.+)/ ) {
				my $header = $1;
				if ( $header =~ /\.(\d+?)\.fasta/ ) {
					$line           = ">$1\n";
					$printing       = defined( $markertaxa{$1} ) ? 0 : 1;
					$markertaxa{$1} = 1;
				}
			}
			print OUTALN $line if $printing;
		}
		close INALN;
		close OUTALN;
		`mv $alignment.fixed $alignment`;
	}

	sub create_temp_read_fasta {
		my %args = @_;
		my $file = $args{file} // miss("file");
		open( TMPREAD, ">$file.tmpread.fasta" );
		print TMPREAD ">blahblahblah\n";
		open( ALNIN, "$file.fasta" );
		my $line = <ALNIN>;
		while ( $line = <ALNIN> ) {
			last if $line =~ /^>/;
			print TMPREAD $line;
		}
	}

	sub filter_fasta {
		my %args         = @_;
		my $input_fasta  = $args{input_fasta} // miss("input_fasta");
		my $output_fasta = $args{output_fasta} // miss("output_fasta");
		my $keep_taxa    = $args{keep_taxa} // miss("keep_taxa");

		# create a pruned fasta
		open( FASTA,       $input_fasta )        || croak "Unable to read $input_fasta";
		open( PRUNEDFASTA, ">" . $output_fasta ) || croak "Unable to write to $output_fasta";
		my $printing = 0;
		while ( my $line = <FASTA> ) {
			chomp $line;
			if ( $line =~ /^>(.+)/ ) {
				$printing = defined( $keep_taxa->{$1} ) ? 1 : 0;
			}
			print PRUNEDFASTA "$line\n" if $printing;
		}
	}

	sub prune_marker {
		my %args         = @_;
		my $tre          = $args{tre} // miss("tre");
		my $distance     = $args{distance} // miss("distance");
		my $fasta        = $args{fasta} // miss("fasta");
		my $pruned_fasta = $args{pruned_fasta} // miss("pruned_fasta");

		#	my $prune_cmd = "$Phylosift::Utilities::pda -k 20000 -g -minlen $args{distance} $args{tre} $args{tre}.pruning.log";
		my $prune_cmd = "pda -k 20000 -g -minlen $distance $tre $tre.pruning.log";
		system("$prune_cmd");

		# read the list of taxa to keep
		open( PRUNE, "$tre.pruning.log" );
		my $intaxa = 0;
		my %keep_taxa;
		while ( my $line = <PRUNE> ) {
			$intaxa = 1 if ( $line =~ /optimal PD set has/ );
			next unless $intaxa;
			chomp $line;
			$keep_taxa{$line} = 1;
			last if ( length($line) < 2 );    # taxa set ends with empty line
		}

		# create a pruned alignment fasta
		filter_fasta( input_fasta => $fasta, output_fasta => $pruned_fasta, keep_taxa => \%keep_taxa );
	}

	sub pd_prune_markers() {
		my %args       = @_;
		my $marker_dir = $args{marker_directory} // miss("marker_directory");

		# prune distance is different for AA and DNA
		# these distances are in substitutions per site
		my %PRUNE_DISTANCE = ( 0 => 0.01, 1 => 0.003 );
		my $REPS_DISTANCE = 0.05;    # reps can be further diverged since we care only about similarity search and not read placement
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {    # zero for aa, one for dna
			foreach my $marker (@markerlist) {
				next unless Phylosift::Utilities::is_protein_marker( marker => $marker );
				my $tre = get_fasttree_tre_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
				my $fasta  = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
				my $pruned = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 1 );
				my $reps_fasta = "";
				my $clean_reps;
				my $pruned_reps;
				prune_marker( distance => $PRUNE_DISTANCE{$dna}, tre => $tre, fasta => $fasta, pruned_fasta => $pruned );

				# if on protein, also create a pruned reps file. Prune more aggressively for this.
				if ( $dna == 0 ) {
					$clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
					$pruned_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1, pruned => 1 );
					prune_marker( distance => $REPS_DISTANCE, tre => $tre, fasta => $clean_reps, pruned_fasta => $pruned_reps );
				}
			}
		}
	}

	sub reconcile_with_ncbi() {
		my %args        = @_;
		my $self        = $args{self} // miss("self");
		my $results_dir = $args{results_directory} // miss("results_directory");
		my $marker_dir  = $args{marker_directory} // miss("marker_directory");
		my $pruned      = $args{pruned} // miss("pruned");
		my $codon_fasta = get_fasta_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_fasta    = get_fasta_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
		my $codon_tre   = get_fasttree_tre_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_tre      = get_fasttree_tre_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
		my $codon_log   = get_fasttree_log_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
		my $aa_log      = get_fasttree_log_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
		open( RECONCILESCRIPT, ">/tmp/ps_reconcile.sh" );
		print RECONCILESCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

# first do the AA tree
taxit create -a "Aaron Darling" -d "simple package for reconciliation only" -l \$1 -f $aa_fasta -t $aa_tre -s $aa_log -Y FastTree -P \$1.updated
pplacer -c \$1.updated -p \$1.updated.tmpread.fasta
# readconciler uses a pplacer tree from a .jplace file to parse out the branch numbers
mangler.pl < \$1.updated.tmpread.jplace > \$1.updated.tmpread.jplace.mangled
readconciler ncbi_tree.updated.tre \$1.updated.tmpread.jplace.mangled marker_taxon_map.updated.txt \$1.updated.taxonmap
rm \$1.updated.tmpread.jplace \$1.updated.tmpread.jplace.mangled \$1.updated.tmpread.fasta \$1.updated.fasttree.log

# then do the codon tree
taxit create -a "Aaron Darling" -d "simple package for reconciliation only" -l \$1 -f $codon_fasta -t $codon_tre -s $codon_log -Y FastTree -P \$1.codon.updated
pplacer -c \$1.codon.updated -p \$1.codon.updated.tmpread.fasta
mangler.pl < \$1.codon.updated.tmpread.jplace > \$1.codon.updated.tmpread.jplace.mangled
readconciler ncbi_tree.updated.tre \$1.codon.updated.tmpread.jplace.mangled marker_taxon_map.updated.txt \$1.codon.updated.taxonmap
rm \$1.codon.updated.tmpread.jplace \$1.codon.updated.tmpread.jplace.mangled \$1.codon.updated.tmpread.fasta \$1.codon.updated.fasttree.log
EOF
		`chmod 755 /tmp/ps_reconcile.sh`;
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		my @jobids;

		foreach my $marker (@markerlist) {

			# create some read files for pplacer to place so we can get its jplace
			create_temp_read_fasta( file => "$marker.updated" );
			create_temp_read_fasta( file => "$marker.codon.updated" );

			# run reconciliation on them
			my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/ps_reconcile.sh $marker";
			my $job      = `$qsub_cmd`;
			$job =~ /Your job (\d+) /;
			push( @jobids, $1 );
		}
		wait_for_jobs( job_ids => @jobids );
		`rm ps_reconcile.sh.o*`;
		`rm ps_reconcile.sh.e*`;
	}

=head2 update_rna

RNA markers contain sequences for many taxa that aren't available in genome databases.
We want the updated markers to retain these original sequences.

=cut

	sub update_rna {
		my %args        = @_;
		my $self        = $args{self} // miss("self");
		my $marker_dir  = $args{marker_dir} // miss("marker_dir");
		my @marker_list = Phylosift::Utilities::gather_markers();
		foreach my $marker (@marker_list) {
			next if Phylosift::Utilities::is_protein_marker( marker => $marker );
			my $base_alignment = Phylosift::Utilities::get_marker_aln_file( $self, $marker );

			#		my $base_alignment = Phylosift::Utilities::get_marker_aln_file(self=>$self, marker=>$marker);
			my $updated_alignment = "$marker.updated.fasta";

			# need to add to gene_ids_aa.txt
			# use greengenes id for taxon ID
			my $aaid_file = "$marker_dir/gene_ids.aa.txt";
			my $max_id    = -1;

			#
			# get the largest gene ID
			open( AAID, $aaid_file ) || croak "Unable to read file $aaid_file";
			while ( my $line = <AAID> ) {
				chomp $line;
				my ( $marker, $tid, $gid ) = split( /\t/, $line );
				$max_id = $gid if $gid > $max_id;
			}
			close AAID;

			#
			# now start assigning IDs to base marker sequences and adding them to the update
			$max_id++;
			open( AAID,   ">>" . $aaid_file )         || croak "Unable to write to $aaid_file";
			open( INALN,  $base_alignment )           || croak "Unable to read $base_alignment";
			open( OUTALN, ">>" . $updated_alignment ) || croak "Unable to write $updated_alignment";
			while ( my $line = <INALN> ) {
				chomp $line;
				if ( $line =~ /^>(.+)/ ) {
					my $countstring = sprintf( "%010u", $max_id );
					my $ggid = $1;
					$ggid =~ s/\/\d+-\d+//g;    # remove trailing bioperl rubbish
					print AAID "$marker\tgg_$ggid\t$countstring\n";
					$line = ">$countstring";
					$max_id++;
				}
				print OUTALN "$line\n";
			}
			close AAID;
		}
	}

=head2 make_codon_submarkers

Find groups of closely related taxa that can would be better analyzed with DNA sequence rather than protein

=cut

	sub make_codon_submarkers {
		my %args = @_;
		my $marker_dir = $args{marker_dir} // miss("marker_dir");

		# this is the maximum distance in amino acid substitutions per site that are allowed
		# on any branch of the tree relating members of a group
		# TODO: tune this value
		my $max_aa_branch_distance = 0.2;
		my ( %id_to_taxon, %marker_taxon_to_id ) = read_gene_ids( file => "$marker_dir/gene_ids_aa.txt" );
		my @subalignments;    # list of subalignments created that will later need marker packages made

		# the following file will provide a mapping of AA gene ID to submarker
		open( SUBTABLE, ">submarkers.txt" ) || croak "Unable to create submarker table file";
		my @marker_list = Phylosift::Utilities::gather_markers();
		unshift( @marker_list, "concat" );
		foreach my $marker (@marker_list) {
			next unless Phylosift::Utilities::is_protein_marker( marker => $marker );
			my $aa_tree = get_fasttree_tre_filename( marker => $marker, dna => 0, updated => 1, pruned => 0 );
			my $codon_alignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0 );
			my $alnio = Phylosift::Utilities::open_SeqIO_object( file => $codon_alignment );
			my $aln = $alnio->next_aln();

			# get groups of taxa that are close
			my @group_table = `segment_tree $aa_tree $max_aa_branch_distance`;
			my %gene_groups;
			my $group_id = 0;
			foreach my $group_line (@group_table) {
				chomp($group_line);
				my @gene_ids = split( /\t/, $group_line );

				# create a new subalignment file
				my $subalignment = $codon_alignment . ".sub$group_id";
				open( SUBALN, ">$subalignment" ) || croak "Unable to write $subalignment";

				# map the gene ID from aa tree into the corresponding ID in the codon data
				foreach my $gene (@gene_ids) {
					my $taxon     = $id_to_taxon{$gene};
					my @codon_ids = $marker_taxon_to_id{$marker}{$taxon};

					# write each sequence into the subalignment
					foreach my $id (@codon_ids) {
						foreach my $seq ( $aln->each_seq_with_id($id) ) {
							print SUBALN ">$id\n";
							print SUBALN $seq->seq() . "\n";
						}
					}

					# add the mapping from gene ID to submarker to the table
					print SUBTABLE "$gene\t$marker\t$group_id\n";
				}
				close SUBALN;
				$group_id++;
			}
		}
		close SUBTABLE;
	}

=head2 make_constrained_tree

Makes a tree for a marker that respects topological constraints given by another marker's tree

=cut

	sub make_constrained_tree {
		my %args              = @_;
		my $constraint_marker = $args{constraint_marker} // miss("constraint_marker");
		my $target_marker     = $args{target_marker} // miss("target_marker");
		my $marker_dir        = $args{marker_dir} // miss("marker_dir");
		my ( %id_to_taxon, %marker_taxon_to_id ) = read_gene_ids("$marker_dir/gene_ids_aa.txt");
		my $constraint_tree = get_fasttree_tre_filename( marker => $constraint_marker, dna => 0, updated => 1, pruned => 1 );

		# first convert the concat tree leaf names to match the 16s
		my $tree = Bio::Phylo::IO->parse(
										  '-file'   => $constraint_tree,
										  '-format' => 'newick',
		)->first;
		my @missing_taxa;    # gene IDs for taxa missing in the 16s
		foreach my $node ( @{ $tree->get_entities } ) {

			# skip this one if it is not a leaf
			next if ( scalar( $node->get_children() ) > 0 );
			my $name = $node->get_name;
			croak "Unable to find taxon for gene $name in marker $constraint_marker" unless defined $id_to_taxon{$name};
			my $taxon_id = $id_to_taxon{$name};

			# now find the taxon in the target data, if it exists
			my $rename;
			if ( defined( $marker_taxon_to_id{$target_marker}{$taxon_id} ) ) {
				my @tids = $marker_taxon_to_id{$target_marker}{$taxon_id};
				$rename = $tids[0];
			} else {
				push( @missing_taxa, $name );
			}
			$node->set_name($rename);
		}
		my $renamed_tree = "$constraint_tree.renamed";
		open( RENAMETREE, ">$renamed_tree" ) || croak "Unable to write to $renamed_tree";
		print RENAMETREE unparse( '-phylo' => $tree, '-format' => 'newick' ) . "\n";
		close RENAMETREE;

		#
		# enumerate the splits in the tree
		my $constraint_splits = "$renamed_tree.splits";
		my $constraint_cl     = "printsplits $renamed_tree $constraint_splits";
		system($constraint_cl);

		#
		# create an amended alignment with any taxa missing from the target data
		my $target_alignment         = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => 1 );
		my $target_alignment_amended = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => 1 ) . ".amended";
		my $column_count             = 0;
		`cp $target_alignment $target_alignment_amended`;
		open( AMENDALN, ">>$target_alignment_amended" ) || croak "Unable to append to $target_alignment_amended\n";
		foreach my $missing (@missing_taxa) {
			print AMENDALN ">$missing\n" . ( "-" x $column_count ) . "\n";
		}
		close AMENDALN;

		#
		# now make a target tree that respects protein tree constraints
		my $target_constrained_tree = get_fasttree_tre_filename( marker => $target_marker, dna => 0, updated => 1, pruned => 1 ) . ".constrained";
		my $target_constrained_log = get_fasttree_log_filename( marker => $target_marker, dna => 0, updated => 1, pruned => 1 ) . ".constrained";
		my $tree_build_cl =
"$Phylosift::Utilities::fasttree -nt -gtr -constraints $constraint_splits -constraintWeight 100 -log $target_constrained_log  $target_alignment > $target_constrained_tree";
		system($tree_build_cl);

		# finally, make a pplacer package with the constrained tree
		my $taxit_cl =
"taxit create -a \"Aaron Darling\" -d \"topology-constrained marker $target_marker\" -l $target_marker -f $target_alignment_amended -t $target_constrained_tree -s $target_constrained_log -Y FastTree -P $target_marker.constrained";
		system($taxit_cl);
	}

=head2 join_trees

Create protein & rna trees with compatible topologies

=cut

	sub join_trees {
		my %args       = @_;
		my $marker_dir = $args{marker_dir} // miss("marker_dir");

		# first apply protein constraints to 16s
		# then apply 16s constraints to protein
		make_constrained_tree( constraint_marker => "concat",       target_marker => "16s_bac_reps", marker_dir => $marker_dir );
		make_constrained_tree( constraint_marker => "16s_bac_reps", target_marker => "concat",       marker_dir => $marker_dir );

		# with a little luck we can now
	}

	sub get_gene_id_file {
		my %args = @_;
		my $dna  = $args{dna} // miss("dna");
		return "gene_ids.codon.txt" if $dna;
		return "gene_ids.aa.txt";
	}

	sub add_taxit_taxonomy {
		my %args       = @_;
		my $marker_dir = $args{marker_dir} // miss("marker_dir");
		my $pruned     = $args{pruned} || 1;

		# only need to run this if it hasnt been done already
		# -- or do we need to re-do it for every taxonomy update?
		#	my $newdb_cl = "taxit new_database -d taxonomy.db";
		#	system($newdb_cl);
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		foreach my $marker (@markerlist) {
			for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {

				# create a taxon id list for this marker
				my $gene_id_file = get_gene_id_file( dna => $dna );
				open( AAIDS,   $gene_id_file )   || croak "Unable to read $gene_id_file";
				open( TAXIDS,  ">tax_ids.txt" )  || croak "Unable to write tax_ids.txt";
				open( SEQINFO, ">seq_info.csv" ) || croak "Unable to write seq_info.csv";
				while ( my $line = <AAIDS> ) {
					chomp $line;
					my @dat = split( /\t/, $line );
					next unless $dat[0] eq $marker;
					print TAXIDS $dat[1] . "\n";
					print SEQINFO "$dat[2],$dat[0],$dat[1],\n";
				}
				close TAXIDS;
				close SEQINFO;
				my $taxtable_cl = "taxit taxtable -d taxonomy.db -t tax_ids.txt -o taxa.csv";
				system($taxtable_cl);

				# gather filenames to stuff into marker package
				my $fasta = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned );
				my $tre = get_fasttree_tre_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned );
				my $log = get_fasttree_log_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned );
				my $taxit_cl =
"taxit create -l $marker.updated -P $marker.updated --taxonomy taxa.csv --seq-info seq_info.csv --tree-stats $log --tree-file $tre --aln-fasta $fasta ";
				system($taxit_cl);
			}
		}
	}

	sub package_markers($) {
		my $marker_dir = shift // miss("marker_dir");
		my @markerlist = Phylosift::Utilities::gather_markers();
		unshift( @markerlist, "concat" );
		foreach my $marker (@markerlist) {

			# move in reps
			my $pruned_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1, pruned => 1 );
			`mv $pruned_reps $marker.updated/$marker.reps`;
			for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {    # zero for aa, one for dna
				                                         # move in taxonmap
				my $taxonmap = get_taxonmap_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
				`mv $taxonmap $marker.updated/$marker.taxonmap`;
			}
		}
		chdir( $marker_dir . "/../" );
		my @timerval = localtime();
		my $datestr  = Phylosift::Utilities::get_date_YYYYMMDD;
		`tar czf markers_$datestr.tgz markers`;
		`rm -f markers.tgz`;
		`ln -s markers_$datestr.tgz markers.tgz`;
	}
	1;
