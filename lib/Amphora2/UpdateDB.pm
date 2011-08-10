package Amphora2::UpdateDB;

use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;

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
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/bacteria.details.txt", $newgenomes);
	get_ebi_from_list("http://www.ebi.ac.uk/genomes/archaea.details.txt", $newgenomes);
#	get_ebi_from_list("http://www.ebi.ac.uk/genomes/eukaryota.details.txt");
}

sub get_ebi_from_list(){
	my $list_url = shift;	# URL to EBI's table of genome characteristics
	my $newgenomes = shift;	# array reference to append file paths of newly downloaded genomes
	`wget $list_url -O list.txt`;


	open( DETAILS, "list.txt" );
	my $line = <DETAILS>;
	while( $line = <DETAILS> ){
		chomp $line;
		my ($acc, $version, $date, $taxid, $description) = split( /\t/, $line );
		my $outfile = $acc;
		$outfile =~ s/\./_/g;
		$outfile .= ".$taxid";
		# if it exists and is newer, don't re-download
		if(-e $outfile){
			my $mtime = (stat($outfile))[9];
			my @timerval=localtime(time);
			my $datestr = (1900+$timerval[5]).($timerval[4]+1).$timerval[3];
			next unless $datestr < $date;
		}
		# either we don't have this one yet, or our version is out of date
		open( WGETTER, "| wget -i - -O $outfile.embl ");
		print WGETTER "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
		close WGETTER;
		my $seq_in = Bio::SeqIO->new(-file => "$outfile.embl");
	
		my $seq_out = Bio::SeqIO->new('-file' => ">$outfile.fasta",
		                               '-format' => "fasta");
		while (my $inseq = $seq_in->next_seq) {
			$seq_out->write_seq($inseq);
		}
		push(@{$newgenomes}, "$outfile.fasta" );
		`rm $outfile.embl`;
	}
	`rm list.txt`;
}

sub get_ncbi_draft_genomes(){
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=fna.tgz ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	my $draft_db = "/share/eisen-d2/amphora2/";
	chdir($draft_db);
	`$ncbi_wget_cmd`;
}

sub find_new_genomes($$$){
	my $genome_dir = shift;
	my $results_dir = shift;
	my $files = shift;
	open( FINDER, "find $genome_dir |" );
	while( my $genome = <FINDER> ){
		chomp $genome;
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
		my $jobid = `qsub -q all.q -q eisen.q /home/koadman/bin/a2sge.sh $file`;
		$jobid =~ /Your job (\d+) /;
		push(@jobids, $1);
	}

	# wait for all jobs to complete
	foreach my $jobid(@jobids){
		while(1){
			my $output = `qstat -j $jobid`;
			last if $output =~ /Following jobs do not exist/;
			sleep(20);
		}
	}
}

sub collate_markers($$$){
	my $results_dir = shift;
	my $marker_dir = shift;
	my $markers = shift;
	foreach my $marker(@{$markers}){
		`rm $marker_dir/$marker.aln`;
		open( FINDER, "find $results_dir -name $marker.aln_hmmer3.trim |" );
		while( my $file = <FINDER> ){
			chomp $file;
			`cat $file >> $marker_dir/$marker.aln`;
		}
		`rm $marker_dir/$marker.codon.aln`;
		open( FINDER, "find $results_dir -name $marker.aln_hmmer3.trim.ffn |" );
		while( my $file = <FINDER> ){
			chomp $file;
			`cat $file >> $marker_dir/$marker.codon.aln`;
		}
	}
}

sub build_marker_trees($){
	# TODO: run raxml on these
	my $dna_tree = "raxmlHPC -m GTRGAMMA -n test -s nucleotides.phy";
	my $aa_tree = "raxmlHPC -m PROTGAMMAWAGF -n test -s amino_acids.phy";
}

sub reconcile_with_ncbi($){
	# TODO: use readconciler
}
