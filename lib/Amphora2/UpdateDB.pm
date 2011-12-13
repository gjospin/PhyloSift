package Amphora2::UpdateDB;

use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Carp;
use Amphora2::Summarize;

=head1 NAME

Amphora2::UpdateDB - Functionality to download new genomes and update the marker database

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


sub get_ebi_genomes($$){
	my $directory = shift;
	my $newgenomes = shift;
	chdir($directory);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/organelle.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/virus.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/phage.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/archaea.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/archaealvirus.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/bacteria.details.txt", $newgenomes);
#	get_ebi_from_list("http://www.ebi.ac.uk/genomes/eukaryota.details.txt");
}

sub get_ebi_from_list(){
	my $list_url = shift;	# URL to EBI's table of genome characteristics
	`wget $list_url -O list.txt`;

	open( DETAILS, "list.txt" );
	my $line = <DETAILS>;
	while( $line = <DETAILS> ){
		# each line in the list file records a genome and some metadata
		chomp $line;
		my ($acc, $version, $date, $taxid, $description) = split( /\t/, $line );
		my $outfile = $acc;
		$outfile =~ s/\./_/g;
		$outfile .= ".$taxid";
		print STDERR "$outfile.fasta\n";
		# if it exists and is newer, don't re-download
		my $download = 1;
		if(-e "$outfile.embl"){
			my $mtime = (stat("$outfile.embl"))[9];
			my @timerval=localtime($mtime);
			my $datestr = (1900+$timerval[5]);
			$datestr .= 0 if $timerval[4] < 9; 
			$datestr .= ($timerval[4]+1);
			$datestr .= 0 if $timerval[3] < 9; 
			$datestr .= $timerval[3];
			next unless int($datestr) < int($date);
			`rm $outfile.fasta`;
		}
		# either we don't have this one yet, or our version is out of date
		open( WGETTER, "| wget -i - -O $outfile.embl ");
		print WGETTER "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
		close WGETTER;
#		eval{
			my $seq_in = Bio::SeqIO->new(-file => "$outfile.embl");	
			my $seq_out = Bio::SeqIO->new('-file' => ">$outfile.fasta",
				                       '-format' => "fasta");
			while (my $inseq = $seq_in->next_seq) {
				$seq_out->write_seq($inseq);
			}
#			`rm $outfile.embl`;
#		} or do {
#			carp "Error processing $outfile.embl\n";
#		}
	}
	`rm list.txt`;
}

sub get_taxid_from_gbk($){
	# first get the taxon ID
	my $file = shift;
	open(GBK, $file);
	my $taxid;
	while( my $l2 = <GBK> ){
		if($l2 =~ /\/db_xref\=\"taxon\:(\d+)/){
			$taxid = $1;
			last;
		}
	}
	return $taxid;
}

sub get_ncbi_finished_genomes($){
	my $directory = shift;
	`mkdir -p $directory`;
	# First download all finished bacterial genomes
	# then for each genome, concatenate all replicons into a FastA file
	# Name the FastA file with the organism and its taxon ID.
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria";
#	`$ncbi_wget_cmd`;
	open( FINDER, "find ftp.ncbi.nih.gov/genomes/Bacteria -type d |" );
	while( my $line = <FINDER> ){
		chomp $line;
		$line =~ /\/Bacteria\/(.+)/;
		my $orgname = $1;
		next if( $orgname =~ /\// );	# could be a malformed genbank directory
		next unless length($1)>1;
		my $seq_out;
		my $fasta_name;
		open( LSSER, "ls $line/*gbk |" );
		while( my $gbk = <LSSER> ){
			chomp $gbk;
			if(!defined($fasta_name)){
				my $taxid = get_taxid_from_gbk($gbk);
				$fasta_name = "$orgname.$taxid.fasta";
				# skip this one if it exists and hasn't been updated
				# POSSIBLE BUG: if only part of this organism's genbank record is updated
				if( -e $fasta_name && (stat($fasta_name))[9] > (stat($gbk))[9]){
					print STDERR "Already have $fasta_name\n";
					last;
				}
				
				$seq_out = Bio::SeqIO->new('-file' => ">$fasta_name",'-format' => "fasta");
			}
			my $seq_in = Bio::SeqIO->new(-file => "$gbk");	
			while (my $inseq = $seq_in->next_seq) {
				$seq_out->write_seq($inseq);
			}
		}
	}
}

sub get_ncbi_draft_genomes($){
	my $directory = shift;
	`mkdir -p $directory`;
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	$ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=\"*fna.tgz\",\"*fna.[0-9].tgz\" ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	open( FINDER, "find . -name \"*.gbk\" |" );
	while( my $line = <FINDER> ){
		chomp $line;
		my $taxid = get_taxid_from_gbk($line);
		# then unpack all the files and combine them
		my $fasta_out = basename($line, ".gbk").".$taxid.fasta";
		my $fna = $line;
		$fna =~ s/\.gbk/\.scaffold\.fna\.tgz/g;
		$fna =~ s/\.scaffold\.fna\.tgz/\.contig\.fna\.tgz/g unless (-e $fna);
		$fna =~ s/\.contig\.fna\.tgz/\.contig\.fna\.1\.tgz/g unless (-e $fna);
		unless( -e $fna ){
			warn "Missing FastA data for $line\n";
			next;
		}elsif($fna =~ /contig\.fna\.1\.tgz/){
			# if the contigs are broken into many files, unpack them all at once
			$fna =~ s/\.contig\.fna\.1\.tgz/\.contig\.fna\.\*\.tgz/g;
		}
		# skip this one if it exists and hasn't been updated
		if( -e $fasta_out && (stat($fasta_out))[9] > (stat($line))[9]){
			print STDERR "Already have $fasta_out\n";
			next;
		}
		# unpack the nucleotide tarball and cat the scaffolds/contigs.
		my @tarfiles = `tar xvzf $fna`;
		my $catline = join(" ",@tarfiles);
		$catline =~ s/\n//g;
		`rm -f $fasta_out`;
		`cat $catline >> $fasta_out`;
		`rm $catline`;
	}
}

sub find_new_genomes($$$){
	my $genome_dir = shift;
	my $results_dir = shift;
	my $files = shift;
	open( FINDER, "find $genome_dir |" );
	while( my $genome = <FINDER> ){
		chomp $genome;
		next unless $genome =~ /\.fasta/;
		my $gbase = basename $genome;
		if( -e "$results_dir/$gbase/alignDir/concat.fasta" ){
			my $ctime = (stat("$results_dir/$gbase/alignDir/concat.fasta"))[9];
			my $mtime = (stat($genome))[9];
			push(@{$files}, $genome) if($ctime < $mtime);
		}else{
			push(@{$files}, $genome);
		}
	}
}

sub qsub_updates($$){
	my $results_dir = shift;	
	my $files = shift;
	my @jobids;
	`mkdir -p $results_dir`;
	foreach my $file(@{$files}){
		chdir($results_dir);
		my $job = `qsub -q all.q -q eisen.q /home/koadman/bin/a2sge.sh $file`;
		$job =~ /Your job (\d+) /;
		push(@jobids, $1);
	}
	wait_for_jobs(@jobids);
}

sub wait_for_jobs{
	# wait for all jobs to complete
	while(my $jobid = shift){
		while(1){
			my $output = `qstat -j $jobid 2>&1`;
			last if $output =~ /Following jobs do not exist/;
			sleep(20);
		}
	}
}

sub collate_markers($$){
	my $results_dir = shift;
	my $marker_dir = shift;
	# get list of markers
	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);
	print STDERR "Markers are ".join(" ",@markerlist)."\n";
	# get a list of genomes available in the results directory
	# this hopefully means we touch each inode over NFS only once
	# NFS is slooooow...
	print STDERR "Working with ".scalar(@markerlist)." markers\n";
	print STDERR "Listing all files in results dir\n";
	my @alldata = `find $results_dir -name "*.aln_hmmer3.trim"`;
	print STDERR "Found ".scalar(@alldata)." files\n";
	foreach my $marker(@markerlist){
		my $cat_ch = ">";
		# find all alignments with this marker
		my @catfiles=();
		foreach my $file(@alldata){
			next unless $file =~ /(.+\.fasta)\/alignDir\/$marker.aln_hmmer3.trim/;
			my $genome = $1;
			chomp($file);
			push(@catfiles, $file);
			# cat the files into the alignment in batches
			# this cuts down on process spawning and file I/O
			# can't do more than about 4k at once because argument lists get too long
			if(@catfiles > 200){
				my $catline = join(" ", @catfiles);
				print STDERR "Found 200 files for marker $marker\n";
				`cat $catline $cat_ch $marker_dir/$marker.updated.fasta`;
				$catline =~ s/\.aln_hmmer3\.trim/\.aln_hmmer3\.trim\.ffn/g;
				`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
				@catfiles = ();
				$cat_ch = ">>";	
			}
		}
		if(@catfiles > 0){
			my $catline = join(" ", @catfiles);
			print STDERR "Last cat for marker $marker\n";
			`cat $catline $cat_ch $marker_dir/$marker.updated.fasta`;
			$catline =~ s/\.aln_hmmer3\.trim/\.aln_hmmer3\.trim\.ffn/g;
			`cat $catline $cat_ch $marker_dir/$marker.codon.updated.fasta`;
		}
		fix_names_in_alignment("$marker_dir/$marker.updated.fasta");
		fix_names_in_alignment("$marker_dir/$marker.codon.updated.fasta");
	}
	`cat $results_dir/*/alignDir/concat.fa >> concat.fa`;
}


# assign 10 digit unique identifiers to each gene
# todo: will need to move into alphabetical space here
sub assign_seqids($$){
	my $results_dir = shift;
	my $marker_dir = shift;
	# get list of markers
	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);
	my $aa_counter = 0;
	my $codon_counter = 0;
	open(my $aa_idtable, ">gene_ids.aa.txt");
	open(my $codon_idtable, ">gene_ids.codon.txt");
	foreach my $marker(@markerlist){
		my $alignment = "$marker.updated.fasta";
		$aa_counter = assign_seqids_for_marker($marker, "$marker.updated.fasta", $aa_idtable, $aa_counter);
		$codon_counter = assign_seqids_for_marker($marker, "$marker.codon.updated.fasta", $codon_idtable, $codon_counter);
	}
}

sub assign_seqids_for_marker($$$){
	my $marker = shift;
	my $alignment = shift;
	my $idtable = shift;
	my $counter = shift;
	open( INALN, $alignment );
	open( OUTALN, ">$alignment.seqids" );
	while( my $line = <INALN> ){
		chomp $line;
		if($line =~ /^>(.+)/){
			my $header = $1;
			my $countstring = sprintf("%010u", $counter);
			print $idtable "$marker\t$header\t$countstring\n";
			$line = ">$countstring";
			$counter++;
		}
		print OUTALN $line."\n";
	}
	close INALN;
	close OUTALN;
	`mv $alignment $alignment.backup`;
	`mv $alignment.seqids $alignment`;
	return $counter;
}

sub get_marker_list{
	my $marker_dir = shift;
	my @markerlist = `find $marker_dir -name \"*.hmm\"`;
	for(my $i=0; $i<@markerlist; $i++){
		chomp $markerlist[$i];
		$markerlist[$i]=~s/\.hmm//g;
		$markerlist[$i]=substr($markerlist[$i], length($marker_dir));
		$markerlist[$i]=~s/^\///g;
	}
	return @markerlist;
}

sub build_marker_trees_fasttree($){
	my $marker_dir = shift;

	open(TREESCRIPT, ">/tmp/a2_tree.sh");
	print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -nt -gtr -log \$1.codon.updated.fasttree.log  \$1.codon.updated.fasta > \$1.codon.updated.tre
/home/koadman/bin/FastTree -log \$1.updated.fasttree.log \$1.updated.fasta > \$1.updated.tre
};
	`chmod 755 /tmp/a2_tree.sh`;

	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);

	my @jobids;
	foreach my $marker(@markerlist){
		next unless (-e "$marker.updated.fasta");
		# run fasttree on them
		my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/a2_tree.sh $marker";
		my $job = `$qsub_cmd`;
		$job =~ /Your job (\d+) /;
		push(@jobids, $1);
	}
	wait_for_jobs(@jobids);
}

sub build_marker_trees_raxml($){
	my $marker_dir = shift;

	open(TREESCRIPT, ">/tmp/a2_tree.sh");
	print TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -pe threaded 3
rm RAxML*\$1*
raxmlHPC -m GTRGAMMA -n \$1.codon.updated -s \$1.codon.updated.phy -T 4
raxmlHPC -m PROTGAMMAWAGF -n \$1.updated -s \$1.updated.phy -T 4
mv RAxML_bestTree.\$1.codon.updated \$1.codon.updated.tre
mv RAxML_bestTree.\$1.updated \$1.updated.tre
mv RAxML_info.\$1.codon.updated \$1.codon.updated.RAxML_info
mv RAxML_info.\$1.updated \$1.updated.RAxML_info
rm RAxML*.\$1*
};
	`chmod 755 /tmp/a2_tree.sh`;

	chdir($marker_dir);
	my @markerlist = get_marker_list($marker_dir);

	my @jobids;
	foreach my $marker(@markerlist){
		# convert to phylip
		next unless (-e "$marker.updated.fasta");
		my $in = Bio::AlignIO->new(-file   => "$marker.updated.fasta", -format => "fasta" );
		my $out = Bio::AlignIO->new(-file   => ">$marker.updated.phy", -format => "phylip" );
		while (my $inseq = $in->next_aln) {
			$out->write_aln($inseq);
		}
		if(-e "$marker.codon.updated.fasta"){
			fix_names_in_alignment("$marker_dir/$marker.codon.updated.fasta");
			my $in_dna = Bio::AlignIO->new(-file   => "$marker.codon.updated.fasta", -format => "fasta" );
			my $out_dna = Bio::AlignIO->new(-file   => ">$marker.codon.updated.phy", -format => "phylip" );
			while (my $inseq = $in_dna->next_aln) {
				$out_dna->write_aln($inseq);
			}
		}
		# run raxml on them
		my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/a2_tree.sh $marker";
		my $job = `$qsub_cmd`;
		$job =~ /Your job (\d+) /;
		push(@jobids, $1);
	}
	wait_for_jobs(@jobids);
}

sub fix_names_in_alignment($){
	my $alignment = shift;
	open( INALN, $alignment );
	open( OUTALN, ">$alignment.fixed" );
	while( my $line = <INALN> ){
		if($line =~ /^>(.+)/){
			my $header = $1;
			if($header =~ /_(\d+)_fasta/){
				$line = ">$1\n";
			}
		}
		print OUTALN $line;
	}
	close INALN;
	close OUTALN;
	`mv $alignment.fixed $alignment`;
}

sub reconcile_with_ncbi($$$){
	my $self = shift;
	my $results_dir = shift;
	my $marker_dir = shift;
	print STDERR "Updating NCBI tree and taxon map...";
#	Amphora2::Summarize::makeNcbiTreeFromUpdate($self, $results_dir, $marker_dir);
	print STDERR "done\n";

	open(RECONCILESCRIPT, ">/tmp/a2_reconcile.sh");
	print RECONCILESCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V

# first do the AA tree
echo ">blahblahblah" > \$1.updated.tmpread.fasta
head -n 2 \$1.updated.fasta | tail -n 1 >> \$1.updated.tmpread.fasta
taxit create -a "Aaron Darling" -d "simple package for reconciliation only" -l \$1 -f \$1.updated.fasta -t \$1.updated.tre -s \$1.updated.fasttree.log -Y FastTree -P \$1
pplacer -c \$1 -p \$1.updated.tmpread.fasta
# readconciler uses a pplacer tree from a .jplace file to parse out the 
perl -p -i -e "s/:(.+?)\\\{(\\d+?)\\\}/\\\{\\\$2\\\}:\\\$1/g" \$1.updated.tmpread.jplace
readconciler ncbi_tree.updated.tre \$1.updated.tmpread.jplace marker_taxon_map.updated.txt \$1.updated.taxonmap

# then do the codon tree
echo ">blahblahblah" > \$1.codon.updated.tmpread.fasta
head -n 2 \$1.codon.updated.fasta | tail -n 1 >> \$1.codon.updated.tmpread.fasta
taxit create -a "Aaron Darling" -d "simple package for reconciliation only" -l \$1 -f \$1.codon.updated.fasta -t \$1.codon.updated.tre -s \$1.codon.updated.fasttree.log -Y FastTree -P \$1
pplacer -c \$1 -p \$1.codon.updated.tmpread.fasta
perl -p -i -e "s/:(.+?)\\\{(\\d+?)\\\}/\\\{\\\$2\\\}:\\\$1/g" \$1.updated.tmpread.jplace
readconciler ncbi_tree.updated.tre \$1.codon.updated.tmpread.jplace marker_taxon_map.updated.txt \$1.codon.updated.taxonmap
EOF

	`chmod 755 /tmp/a2_reconcile.sh`;

	my @markerlist = get_marker_list($marker_dir);
	my @jobids;
	foreach my $marker(@markerlist){
		# run reconciliation on them
		my $qsub_cmd = "qsub -q all.q -q eisen.q /tmp/a2_reconcile.sh $marker";
		my $job = `$qsub_cmd`;
		$job =~ /Your job (\d+) /;
		push(@jobids, $1);
	}
	wait_for_jobs(@jobids);
}

sub package_markers($){
	my $marker_dir = shift;
	chdir($marker_dir."/../");
	
	my @timerval = localtime();
	my $datestr = Amphora2::Utilities::get_date_YYYYMMDD;
	`tar czf markers_$datestr.tgz markers`;
	`rm -f markers.tgz`;
	`ln -s markers_$datestr.tgz markers.tgz`;
}

1;

