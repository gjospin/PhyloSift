package Phylosift::UpdateDB;
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Carp;
use Phylosift::Summarize;

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
	my $directory = shift;
	chdir($directory);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/organelle.details.txt");
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/virus.details.txt");
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/phage.details.txt");
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/archaea.details.txt");
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/archaealvirus.details.txt");
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/bacteria.details.txt");

	#	get_ebi_from_list("http://www.ebi.ac.uk/genomes/eukaryota.details.txt");
}

sub get_ebi_from_list() {
	my $list_url = shift;    # URL to EBI's table of genome characteristics
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
		my $seq_in = Phylosift::Utilities::open_SeqIO_object( file => "$outfile.embl" );
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

sub get_taxid_from_gbk($) {

	# first get the taxon ID
	my $file = shift;
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

sub get_ncbi_finished_genomes($) {
	my $directory = shift;
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
				my $taxid = get_taxid_from_gbk($gbk);
				$fasta_name = "$orgname.$taxid.fasta";

				# skip this one if it exists and hasn't been updated
				# POSSIBLE BUG: if only part of this organism's genbank record is updated
				if ( -e $fasta_name && ( stat($fasta_name) )[9] > ( stat($gbk) )[9] ) {
					print STDERR "Already have $fasta_name\n";
					last;
				}
				$seq_out = Bio::SeqIO->new( '-file' => ">$fasta_name", '-format' => "fasta" );
			}
			my $seq_in = Bio::SeqIO->new( -file => "$gbk" );
			while ( my $inseq = $seq_in->next_seq ) {
				$seq_out->write_seq($inseq);
			}
		}
	}
}

sub get_ncbi_draft_genomes($) {
	my $directory = shift;
	`mkdir -p $directory`;
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	$ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=\"*fna.tgz\",\"*fna.[0-9].tgz\" ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	open( FINDER, "find . -name \"*.gbk\" |" );

	while ( my $line = <FINDER> ) {
		chomp $line;
		my $taxid = get_taxid_from_gbk($line);

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

sub find_new_genomes($$$) {
	my $genome_dir  = shift;
	my $results_dir = shift;
	my $files       = shift;
	open( FINDER, "find $genome_dir |" );
	while ( my $genome = <FINDER> ) {
		chomp $genome;
		next unless $genome =~ /\.fasta/;
		my $gbase = basename $genome;
		if ( -e "$results_dir/$gbase/alignDir/concat.fasta" ) {
			my $ctime = ( stat("$results_dir/$gbase/alignDir/concat.fasta") )[9];
			my $mtime = ( stat($genome) )[9];
			push( @{$files}, $genome ) if ( $ctime < $mtime );
		} else {
			push( @{$files}, $genome );
		}
	}
}

sub qsub_updates($$) {
	my $results_dir = shift;
	my $files       = shift;
	my @jobids;
	`mkdir -p $results_dir`;
	foreach my $file ( @{$files} ) {
		chdir($results_dir);
		my $job = `qsub -q all.q -q eisen.q /home/koadman/bin/pssge.sh $file`;
		$job =~ /Your job (\d+) /;
		push( @jobids, $1 );
	}
	wait_for_jobs(@jobids);
}

sub wait_for_jobs {

	# wait for all jobs to complete
	while ( my $jobid = shift ) {
		while (1) {
			my $output = `qstat -j $jobid 2>&1`;
			last if $output =~ /Following jobs do not exist/;
			sleep(20);
		}
	}
}

sub collate_markers($$) {
	my $results_dir = shift;
	my $marker_dir  = shift;

	# get list of markers
	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);
	print STDERR "Markers are " . join( " ", @markerlist ) . "\n";

	# get a list of genomes available in the results directory
	# this hopefully means we touch each inode over NFS only once
	# NFS is slooooow...
	print STDERR "Working with " . scalar(@markerlist) . " markers\n";
	print STDERR "Listing all files in results dir\n";
	my @alldata = `find $results_dir -name "*.aln_hmmer3.trim"`;
	print STDERR "Found " . scalar(@alldata) . " files\n";
	foreach my $marker (@markerlist) {
		my $cat_ch = ">";

		# find all alignments with this marker
		my @catfiles = ();
		foreach my $file (@alldata) {
			next unless $file =~ /(.+\.fasta)\/alignDir\/$marker.aln_hmmer3.trim/;
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
				$catline =~ s/\.aln_hmmer3\.trim/\.aln_hmmer3\.trim\.ffn/g;
				`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
				@catfiles = ();
				$cat_ch   = ">>";
			}
		}
		if ( @catfiles > 0 ) {
			my $catline = join( " ", @catfiles );
			print STDERR "Last cat for marker $marker\n";
			`cat $catline $cat_ch $marker_dir/$marker.updated.fasta`;
			$catline =~ s/\.aln_hmmer3\.trim/\.aln_hmmer3\.trim\.ffn/g;
			`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
		}
		fix_names_in_alignment("$marker_dir/$marker.updated.fasta");
		fix_names_in_alignment("$marker_dir/$marker.codon.updated.fasta");
	}
	`rm -f $marker_dir/concat.updated.fasta`;
	`rm -f $marker_dir/concat.codon.updated.fasta`;
	`cat $results_dir/*/alignDir/concat.fasta >> $marker_dir/concat.updated.fasta`;
	`cat $results_dir/*/alignDir/concat-dna.fasta >> $marker_dir/concat.codon.updated.fasta`;
	fix_names_in_alignment("$marker_dir/concat.updated.fasta");
	fix_names_in_alignment("$marker_dir/concat.codon.updated.fasta");
}

# assign 10 digit unique identifiers to each gene
# todo: will need to move into alphabetical space here
sub assign_seqids($) {
	my $marker_dir = shift;

	# get list of markers
	chdir($marker_dir);
	my @markerlist    = get_marker_list($marker_dir);
	my $aa_counter    = 0;
	my $codon_counter = 0;
	open( my $aa_idtable,    ">gene_ids.aa.txt" );
	open( my $codon_idtable, ">gene_ids.codon.txt" );
	push( @markerlist, "concat" );
	foreach my $marker (@markerlist) {
		$aa_counter    = assign_seqids_for_marker( $marker, "$marker.updated.fasta",       $aa_idtable,    $aa_counter );
		$codon_counter = assign_seqids_for_marker( $marker, "$marker.codon.updated.fasta", $codon_idtable, $codon_counter );
	}
}

sub assign_seqids_for_marker($$$) {
	my $marker    = shift;
	my $alignment = shift;
	my $idtable   = shift;
	my $counter   = shift;
	open( INALN,  $alignment );
	open( OUTALN, ">$alignment.seqids" );
	my %seen_ids;
	my $printing = 0;

	while ( my $line = <INALN> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			my $header = $1;
			if ( defined( $seen_ids{$header} ) ) {
				$printing = 0;
				next;
			}
			$seen_ids{$header} = 1;
			$printing = 1;
			my $countstring = sprintf( "%010u", $counter );
			print $idtable "$marker\t$header\t$countstring\n";
			$line = ">$countstring";
			$counter++;
		}
		print OUTALN $line . "\n" if $printing;
	}
	close INALN;
	close OUTALN;
	`mv $alignment.seqids $alignment`;
	return $counter;
}

sub get_marker_list {
	my $marker_dir = shift;
	my @markerlist = `find $marker_dir -name \"*.hmm\"`;
	for ( my $i = 0 ; $i < @markerlist ; $i++ ) {
		chomp $markerlist[$i];
		$markerlist[$i] =~ s/\.hmm//g;
		$markerlist[$i] = substr( $markerlist[$i], length($marker_dir) );
		$markerlist[$i] =~ s/^\///g;
	}
	return @markerlist;
}

sub get_marker_name_base {
	my %args = @_;
	my $name = $args{marker};
	$name .= ".codon"   if $args{dna};
	$name .= ".updated" if $args{updated};
	$name .= ".pruned"  if $args{pruned};
	return $name;
}

sub get_fasta_filename {
	my %args = @_;
	my $name = get_marker_name_base(%args);
	$name .= ".fasta";
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

sub build_marker_trees_fasttree($$) {
	my $marker_dir  = shift;
	my $pruned      = shift;
	my $codon_fasta = get_fasta_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_fasta    = get_fasta_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_tre   = get_fasttree_tre_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_tre      = get_fasttree_tre_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_log   = get_fasttree_log_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_log      = get_fasttree_log_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
	open( TREESCRIPT, ">/tmp/ps_tree.sh" );
	print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -nt -gtr -log $codon_log  $codon_fasta > $codon_tre
/home/koadman/bin/FastTree -log $aa_log $aa_fasta > $aa_tre
};
	`chmod 755 /tmp/ps_tree.sh`;
	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);
	unshift( @markerlist, "concat" );
	my @jobids;

	foreach my $marker (@markerlist) {
		next unless ( -e $aa_fasta );

		# run fasttree on them
		my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/ps_tree.sh $marker";
		my $job      = `$qsub_cmd`;
		$job =~ /Your job (\d+) /;
		push( @jobids, $1 );
	}
	wait_for_jobs(@jobids);
	`rm ps_tree.sh.*`;
}

sub build_marker_trees_raxml($) {
	my $marker_dir = shift;
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
	my @markerlist = get_marker_list($marker_dir);
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
			fix_names_in_alignment("$marker_dir/$marker.codon.updated.fasta");
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
	wait_for_jobs(@jobids);
	`rm ps_tree.sh.*`;
}

sub fix_names_in_alignment($) {
	my $alignment = shift;
	my %markertaxa;

	# naive means to remove most trivial duplicates. TODO: allow divergent paralogs to remain.
	my $printing = 1;
	open( INALN,  $alignment );
	open( OUTALN, ">$alignment.fixed" );
	while ( my $line = <INALN> ) {
		if ( $line =~ /^>(.+)/ ) {
			my $header = $1;
			if ( $header =~ /_(\d+)_fasta/ ) {
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

sub createTempReadFasta {
	my $file = shift;
	open( TMPREAD, ">$file.tmpread.fasta" );
	print TMPREAD ">blahblahblah\n";
	open( ALNIN, "$file.fasta" );
	my $line = <ALNIN>;
	while ( $line = <ALNIN> ) {
		last if $line =~ /^>/;
		print TMPREAD $line;
	}
}

sub prune_marker {
	my %args      = @_;
	my $prune_cmd = "pda -k 20000 -g -minlen $args{distance} $args{tre} $args{tre}.pruning.log";
	system("$prune_cmd");

	# read the list of taxa to keep
	open( PRUNE, "$args{tre}.pruning.log" );
	my $intaxa = 0;
	my %keep_taxa;
	while ( my $line = <PRUNE> ) {
		$intaxa = 1 if ( $line =~ /optimal PD set has/ );
		next unless $intaxa;
		chomp $line;
		$keep_taxa{$line} = 1;
		last if ( length($line) < 2 );    # taxa set ends with empty line
	}

	# create a pruned fasta
	open( FASTA,       $args{fasta} );
	open( PRUNEDFASTA, $args{pruned_fasta} );
	my $printing = 0;
	while ( my $line = <FASTA> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			$printing = defined( $keep_taxa{$1} ) ? 1 : 0;
		}
		print PRUNEDFASTA "$line\n" if $printing;
	}
}

sub pd_prune_markers($) {
	my $marker_dir = shift;

	# prune distance is different for AA and DNA
	my %PRUNE_DISTANCE = ( 0 => 0.01, 1 => 0.003 );
	my @markerlist = get_marker_list($marker_dir);
	unshift( @markerlist, "concat" );
	for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {    # zero for aa, one for dna
		foreach my $marker (@markerlist) {
			my $tre = get_fasttree_tre_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
			my $fasta  = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
			my $pruned = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 1 );
			prune_marker( distance => $PRUNE_DISTANCE{$dna}, tre => $tre, fasta => $fasta, pruned_fasta => $pruned );
		}
	}
}

sub reconcile_with_ncbi($$$$) {
	my $self        = shift;
	my $results_dir = shift;
	my $marker_dir  = shift;
	my $pruned      = shift;
	print STDERR "Updating NCBI tree and taxon map...";
	Phylosift::Summarize::makeNcbiTreeFromUpdate( $self, $results_dir, $marker_dir );
	print STDERR "done\n";
	my $codon_fasta = get_fasta_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_fasta    = get_fasta_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_tre = get_fasttree_tre_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_tre    = get_fasttree_tre_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_log = get_fasttree_log_filename( marker => '\\$1', dna => 1, updated => 1, pruned => $pruned );
	my $aa_log    = get_fasttree_log_filename( marker => '\\$1', dna => 0, updated => 1, pruned => $pruned );
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
	my @markerlist = get_marker_list($marker_dir);
	unshift( @markerlist, "concat" );
	my @jobids;

	foreach my $marker (@markerlist) {

		# create some read files for pplacer to place so we can get its jplace
		createTempReadFasta("$marker.updated");
		createTempReadFasta("$marker.codon.updated");

		# run reconciliation on them
		my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/ps_reconcile.sh $marker";
		my $job      = `$qsub_cmd`;
		$job =~ /Your job (\d+) /;
		push( @jobids, $1 );
	}
	wait_for_jobs(@jobids);
	`rm ps_reconcile.sh.o*`;
	`rm ps_reconcile.sh.e*`;
}

sub package_markers($) {
	my $marker_dir = shift;
	chdir( $marker_dir . "/../" );
	my @timerval = localtime();
	my $datestr  = Phylosift::Utilities::get_date_YYYYMMDD;
	`tar czf markers_$datestr.tgz markers`;
	`rm -f markers.tgz`;
	`ln -s markers_$datestr.tgz markers.tgz`;
}
1;
