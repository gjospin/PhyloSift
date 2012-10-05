package Phylosift::UpdateDB;
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use Carp;
use Phylosift::Summarize;
use Phylosift::Settings;
use Phylosift::Utilities;
use Phylosift::MarkerBuild;
use Bio::Phylo::IO qw(parse unparse);
use LWP::Simple;
use FindBin;

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
	my %args = @_;
	my $directory = $args{directory} || miss("directory");
	chdir($directory);
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/organelle.details.txt" );
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/virus.details.txt" );
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/phage.details.txt" );
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaea.details.txt" );
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaealvirus.details.txt" );
#	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/bacteria.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/eukaryota.details.txt");
}

sub get_ebi_from_list {
	my %args = @_;
	my $list_url = $args{url} || miss("url");    # URL to EBI's table of genome characteristics
	`wget $list_url -O "list.txt"`;
	my $DETAILS = ps_open( "list.txt" );
	my $line = <$DETAILS>;
	my %taxa_details; # {$taxon_id}{$accession} = [$date, $version, $size]
	while ( $line = <$DETAILS> ) {
		# each line in the list file records a genome and some metadata
		chomp $line;
		my ( $acc, $version, $date, $taxid, $description ) = split( /\t/, $line );
		next if $description =~ /itochondri/; # skip mitochondrion
		next if $description =~ /hloroplas/; # skip chloroplast
		$taxa_details{$taxid}{$acc} = {dater=>$date, version=>$version};
	}
	foreach my $taxid(keys(%taxa_details)){
		# add up the accession sizes
		my $total_size = 0;
		my $facc;
		my $date;
		
		foreach my $acc(keys(%{$taxa_details{$taxid}})){
			# if it's too big, don't download
			my $url = "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
			my $ua = LWP::UserAgent->new(max_size=>200);
			my $response = $ua->get($url);
			my $content = $response->decoded_content();
			my $sizer = $1 if $content =~ /(\d+) BP/;
			$total_size += $sizer;
			$facc = $acc unless defined $facc;
			print STDERR "taxid $taxid, acc $acc\n";
			$date = $taxa_details{$taxid}{$acc}->{dater} unless defined $date;
		}
		print STDERR "Total genome size of taxid $taxid is $total_size\n";
		
		next if $total_size > 200000000;  # 100Mbp limit
		
		my $outfile = $facc;
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
#			next unless int($datestr) < int($date);
			`rm "$outfile.fasta"`;
		}

		# either we don't have this one yet, or our version is out of date
		my $catch = ">";
		foreach my $acc(keys(%{$taxa_details{$taxid}})){
			my $WGETTER = ps_open("| wget -i - -O $outfile.embl.tmp " );
			print $WGETTER "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
			close $WGETTER;
			`cat $outfile.embl.tmp $catch $outfile.embl ; rm $outfile.embl.tmp`;
			$catch = ">>";
		}

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
	`rm "list.txt"`;
}

sub get_taxid_from_gbk {
	my %args = @_;
	my $file = $args{file} || miss("file");

	# first get the taxon ID
	my $GBK = ps_open( $file );
	my $taxid;
	while ( my $l2 = <$GBK> ) {
		if ( $l2 =~ /\/db_xref\=\"taxon\:(\d+)/ ) {
			$taxid = $1;
			last;
		}
	}
	return $taxid;
}

sub get_ncbi_finished_genomes {
	my %args = @_;
	my $directory = $args{directory} || miss("directory");
	`mkdir -p "$directory"`;

	# First download all finished bacterial genomes
	# then for each genome, concatenate all replicons into a FastA file
	# Name the FastA file with the organism and its taxon ID.
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria";

	#	`$ncbi_wget_cmd`;
	my $FINDER = ps_open( "find ftp.ncbi.nih.gov/genomes/Bacteria -type d |" );
	while ( my $line = <$FINDER> ) {
		chomp $line;
		$line =~ /\/Bacteria\/(.+)/;
		my $orgname = $1;
		next if ( $orgname =~ /\// );    # could be a malformed genbank directory
		next unless length($1) > 1;
		my $seq_out;
		my $fasta_name;
		my $LSSER = ps_open( "ls $line/*gbk |" );

		while ( my $gbk = <$LSSER> ) {
			chomp $gbk;
			if ( !defined($fasta_name) ) {
				my $taxid = get_taxid_from_gbk( file => $gbk );
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

sub get_ncbi_wgs_genomes {
	my %args = @_;
	my $directory = $args{directory} || miss("directory");
	`mkdir -p "$directory"`;

	# First download all finished bacterial genomes
	# then for each genome, concatenate all replicons into a FastA file
	# Name the FastA file with the organism and its taxon ID.
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --limit-size=25000000 --accept=gbff.gz ftp://ftp.ncbi.nih.gov/genbank/wgs/";

#	`$ncbi_wget_cmd`;
	my $FINDER = ps_open( "find ftp.ncbi.nih.gov/genbank/wgs -name \"*.gz\" |" );
	while ( my $line = <$FINDER> ) {
		next if $line =~ /\.mstr\.gbff\.gz/;
		chomp $line;
		`gunzip $line`;
		$line =~ s/\.gz//g;
		my $orgname = $1 if $line =~ /wgs\.(....)\./;
		my $taxid = get_taxid_from_gbk( file => $line );
		my $fasta_name = "$orgname.$taxid.fasta";
		if ( -e $fasta_name && ( stat($fasta_name) )[9] > ( stat($line) )[9] ) {
			print STDERR "Already have $fasta_name\n";
			last;
		}
		my $seq_out = Phylosift::Utilities::open_SeqIO_object( file => ">$fasta_name", format => "FASTA" );
		my $seq_in = Phylosift::Utilities::open_SeqIO_object( file => "$line", format => "genbank" );
		while ( my $inseq = $seq_in->next_seq ) {
			$seq_out->write_seq($inseq);
		}
	}
}

sub get_ncbi_draft_genomes {
	my %args = @_;
	my $directory = $args{directory} || miss("directory");
	`mkdir -p "$directory"`;
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	$ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=\"*fna.tgz\",\"*fna.[0-9].tgz\" ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT";
	`$ncbi_wget_cmd`;
	my $FINDER = ps_open( "find . -name \"*.gbk\" |" );

	while ( my $line = <$FINDER> ) {
		chomp $line;
		my $taxid = get_taxid_from_gbk( file => $line );

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
		my @tarfiles  = `tar xvzf "$fna"`;
		my $tarstatus = ( $? >> 8 );
		my $catline   = join( " ", @tarfiles );
		$catline =~ s/\n//g;
		`rm -f "$fasta_out"`;
		`cat "$catline" >> "$fasta_out"` if ( $tarstatus == 0 );
		`rm "$catline"`;
	}
}

sub find_new_genomes {
	my %args        = @_;
	my $genome_dir  = $args{genome_directory} || miss("genome_directory");
	my $results_dir = $args{results_directory} || miss("results_directory");
	my $files       = $args{files} || miss("files");
	my $FINDER = ps_open( "find $genome_dir |" );
	while ( my $genome = <$FINDER> ) {
		chomp $genome;
		next unless $genome =~ /\.fasta/;
		my $gbase = basename $genome;
		print STDERR "Looking for $results_dir/$gbase/alignDir/concat.fasta\n";
		if ( -e "$results_dir/$gbase/alignDir/concat.fasta" && -s "$results_dir/$gbase/alignDir/concat.fasta" > 0 ) {
			my $ctime = ( stat("$results_dir/$gbase/alignDir/concat.fasta") )[9];
			my $mtime = ( stat($genome) )[9];
			push( @{$files}, $genome ) if ( $ctime < $mtime );
#			push( @{$files}, $genome );
			print STDERR "Found up-to-date $gbase\n" if ( $ctime >= $mtime );
		} else {
			push( @{$files}, $genome );
		}
	}
}

# limit on the number of jobs that can be queued to SGE at once
use constant MAX_SGE_JOBS => 5000;

sub qsub_updates {
	my %args        = @_;
	my $local_directory = $args{local_directory} || miss ("local_directory");
	my $files       = $args{files} || miss("files");
	my $params;
#	$params      = "--updated_markers=0 --extended";
	# qsub a slew of phylosift jobs and have them report results to a local directory on this host for concatenation
	my $hostname = `hostname`;
	chomp $hostname;
	`mkdir -p "$local_directory"`;
	
	my @jobids;
	my $PHYLOSIFTSCRIPT = ps_open( ">/tmp/pssge.sh" );
	print $PHYLOSIFTSCRIPT qq{#!/bin/sh
#\$ -V
#\$ -S /bin/bash
export PATH=\$PATH:/home/koadman/development/PhyloSift/bin/
export PERL5LIB=\$HOME/lib/perl5:/home/koadman/development/PhyloSift/lib
WORKDIR=/state/partition1/koadman/phylosift/\$JOB_ID
mkdir -p \$WORKDIR
if [ \$? -gt 0 ]; then
    WORKDIR=/data/scratch/koadman/phylosift/\$JOB_ID
	mkdir -p \$WORKDIR
fi
cd \$WORKDIR
phylosift search $params --isolate --besthit --unique \$1
phylosift align $params --isolate --besthit --unique \$1
rm -rf PS_temp/*/treeDir
rm -rf PS_temp/*/blastDir
rm -rf PS_temp/*/isolates.fasta
rm -rf PS_temp/*/alignDir/PMP*.newCandidate

#scp -r PS_temp/* $hostname:$local_directory
cp -r PS_temp/* $local_directory
rm -rf \$WORKDIR
};
	my $job_count = 0;

	foreach my $file ( @{$files} ) {
		$job_count++;
		qsub_job(script=>"/tmp/pssge.sh", job_ids=>\@jobids, script_args=>[$file] );

		# check whether we've hit the limit for queued jobs, and rest if needed
		if ( $job_count == MAX_SGE_JOBS ) {
			wait_for_jobs(job_ids=>\@jobids);
			@jobids    = ();
			$job_count = 0;
		}
	}
	wait_for_jobs( job_ids => \@jobids );
}

sub wait_for_jobs {
	my %args = @_;
	my $jobref = $args{job_ids} || miss("job_ids");
	my @job_ids = @$jobref;

	# wait for all jobs to complete
	foreach my $jobid ( @job_ids ) {
		while (1) {
			my $output = `qstat -j $jobid 2>&1`;
			last unless defined($output);
			last if $output =~ /Following jobs do not exist/;
			sleep(20);
		}
	}
}

sub concat_marker_files {
	my %args        = @_;
	my $marker = $args{marker};
	my $local_directory = $args{local_directory} || miss ("local_directory");
	my $cfref = $args{catfiles};
	my $cat_ch = $args{cat_ch};
	
	my $catline = join( " ", @$cfref );
	$catline .= " ";
	my $fasta = get_fasta_filename(marker=>$marker, updated=>1);
	`cat $catline $cat_ch "$local_directory/$fasta"`;
	$catline =~ s/updated\.fasta /codon\.updated\.fasta /g;
	if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
		my $codon_fasta = get_fasta_filename(marker=>$marker, updated=>1, dna=>1);
		`cat $catline $cat_ch "$local_directory/$codon_fasta"`;
	}
	unless($marker eq "concat"){
		$catline =~ s/\.codon\.updated\.fasta/\.unmasked/g;
		my $reps = get_reps_filename( marker => $marker, updated => 1 );
		`cat $catline $cat_ch "$local_directory/$reps"`;
	}
}

sub collate_markers {
	my %args        = @_;
	my $local_directory = $args{local_directory} || miss ("local_directory");
	my $marker_dir  = $args{marker_dir} || miss("marker_dir");
	my $taxon_knockouts = $args{taxon_knockouts};

	# make the directory if it doesn't yet exist
	`mkdir -p $marker_dir`;

	# get list of markers
	chdir($marker_dir);
	my @markerlist = Phylosift::Utilities::gather_markers();
	print STDERR "Markers are " . join( " ", @markerlist ) . "\n";
	
	my %ko_list;
	if(defined($taxon_knockouts)){
		my $KO = ps_open($taxon_knockouts);
		while(my $line = <$KO>){
			chomp $line;
			$ko_list{$line}=1;
		}
	}
	
	my $bs = "/tmp/ps_build_marker.sh";
	write_marker_build_script(batch_script=>$bs);

	# get a list of genomes available in the results directory
	# this hopefully means we touch each inode over NFS only once
	# NFS is slooooow...
	print STDERR "Working with " . scalar(@markerlist) . " markers\n";
	print STDERR "Listing all files in results dir\n";
	my @alldata = `find "$local_directory" -name "*.fasta"`;
	print STDERR "Found " . scalar(@alldata) . " files\n";
	unshift( @markerlist, "concat" );
	my %alltaxa;
	my @job_ids;
	foreach my $marker (@markerlist) {
		my $cat_ch = ">";    # first time through ensures that existing files get clobbered

		# find all alignments with this marker
		my @catfiles = ();
		foreach my $file (@alldata) {
			next unless $file =~ /(.+\.fasta)\/alignDir\/$marker.updated.fasta/;
			my $genome = $1;
			my $taxon = $1 if $genome =~ /\.(\d+)\.fasta/;
			next if( defined($taxon) && defined($ko_list{$taxon}));
			chomp($file);
			push( @catfiles, $file );

			# cat the files into the alignment in batches
			# this cuts down on process spawning and file I/O
			# can't do more than about 4k at once because argument lists get too long
			my $batch_size = 1000;
			if ( @catfiles > $batch_size ) {
				print STDERR "Found $batch_size files for marker $marker\n";
				concat_marker_files(marker=>$marker, local_directory=>$local_directory, catfiles=>\@catfiles, cat_ch=>$cat_ch);
				@catfiles = ();
				$cat_ch   = ">>";    # now add to existing files
			}
		}
		if ( @catfiles > 0 ) {
			print STDERR "Last cat for marker $marker\n";
			concat_marker_files(marker=>$marker, local_directory=>$local_directory, catfiles=>\@catfiles, cat_ch=>$cat_ch);
		}
		
		# now rename sequences with their taxon IDs 
		my $fasta = get_fasta_filename(marker=>$marker, updated=>1);
		next unless -e "$local_directory/$fasta";
		create_taxon_id_table(alignment => "$local_directory/$fasta", output => "$marker_dir/$fasta.taxon_ids", alltaxa=>\%alltaxa);
#		fix_names_in_alignment( alignment => "$local_directory/$fasta" );
		`cp "$local_directory/$fasta" "$marker_dir/$fasta"`;
		my $reps = get_reps_filename( marker => $marker, updated => 1 );
		my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
		if(-e "$local_directory/$reps" ){
			clean_representatives(infile=>"$local_directory/$reps", outfile=>"$local_directory/$clean_reps" );
#			fix_names_in_alignment( alignment => "$local_directory/$clean_reps" );
			`cp "$local_directory/$clean_reps" "$marker_dir/$clean_reps"`;
		}
		debug "Launching marker build for $marker\n";
		my $job_id = launch_marker_build(marker=>$marker, dna=>0, batch_script=>$bs);
		push(@job_ids, $job_id);
		
		my $codon_fasta = get_fasta_filename(marker=>$marker, updated=>1,dna=>1);
		next unless -e "$local_directory/$codon_fasta";
		create_taxon_id_table(alignment => "$local_directory/$codon_fasta", output => "$marker_dir/$codon_fasta.taxon_ids", alltaxa=>\%alltaxa);
##		fix_names_in_alignment( alignment => "$local_directory/$codon_fasta" );
		`cp "$local_directory/$codon_fasta" "$marker_dir/$codon_fasta"`;

		debug "Launching marker build for $marker.codon\b";
		$job_id = launch_marker_build(marker=>$marker, dna=>1, batch_script=>$bs);
		push(@job_ids, $job_id);
	}
	my @taxonids = keys(%alltaxa);
	Phylosift::MarkerBuild::make_ncbi_subtree(out_file=>"$marker_dir/ncbi_tree.updated.tre", taxon_ids=>\@taxonids);

	wait_for_jobs( job_ids => \@job_ids );
}

sub clean_representatives {
	my %args    = @_;
	my $file    = $args{infile} || miss("infile");
	my $outfile = $args{outfile} || miss("outfile");
	my $INALN = ps_open( $file );
	my $OUTALN = ps_open( ">$outfile" );
	while ( my $line = <$INALN> ) {
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
			$line =~ tr/[a-z]/[A-Z]/; # make all chars uppercase for lastal (it uses lowercase as soft-masking)
			$line .= "\n";       # trailing newline was removed above
		}
		print $OUTALN $line;
	}
}

=head2 assign_seqids
 assign 10 digit unique identifiers to each gene
 todo: will need to move into alphabetical space here
=cut

sub assign_seqids {
	my %args = @_;
	my $marker_dir = $args{marker_directory} || miss("marker_directory");

	# get list of markers
	chdir($marker_dir);
	my @markerlist    = Phylosift::Utilities::gather_markers();
	my $aa_counter    = 0;
	my $codon_counter = 0;
	my $AA_IDTABLE = ps_open( ">".get_gene_id_file(dna=>0) );
	my $CODON_IDTABLE = ps_open( ">".get_gene_id_file(dna=>1) );
	push( @markerlist, "concat" );
	foreach my $marker (@markerlist) {
		my %id_mapping       = ();
		my %id_mapping_codon = ();
		my $fasta = get_fasta_filename(marker=>$marker, updated=>1);
		next unless -e $fasta;
		print STDERR "marker $marker seqids\n";
		$aa_counter = assign_seqids_for_marker(
												marker       => $marker,
												alignment    => $fasta,
												IDTABLE      => $AA_IDTABLE,
												counter      => $aa_counter,
												existing_ids => \%id_mapping
		);
		# now put the IDs in the reps file
		print STDERR "marker $marker rep seqids\n";
		unless( $marker =~ /concat/ ){
			my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
			if(-e $clean_reps){
				assign_seqids_for_marker( marker => $marker, alignment => $clean_reps, counter => 0, existing_ids => \%id_mapping );
			}
		}

		my $codon_fasta = get_fasta_filename(marker=>$marker, updated=>1, dna=>1);
		next unless -e $codon_fasta;
		$codon_counter = assign_seqids_for_marker(
												   marker       => $marker,
												   alignment    => $codon_fasta,
												   IDTABLE      => $CODON_IDTABLE,
												   counter      => $codon_counter,
												   existing_ids => \%id_mapping_codon
		);

	}
}

sub assign_seqids_for_marker {
	my %args         = @_;
	my $marker       = $args{marker} || miss("marker");
	my $alignment    = $args{alignment} || miss("alignment");
	my $IDTABLE      = $args{IDTABLE} || 0;
	my $counter      = $args{counter} || 0;
	my $existing_ids = $args{existing_ids} || miss("existing_ids");
	my $INALN = ps_open( $alignment );
	my $OUTALN = ps_open( ">$alignment.seqids" );
	my %mapped_ids;
	my $printing = 0;

	while ( my $line = <$INALN> ) {
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
				print $IDTABLE "$marker\t$header\t$countstring\n" if fileno $IDTABLE;
			}
			$mapped_ids{$header} = $countstring;
			$line = ">$countstring";
			$counter++;
		}
		print $OUTALN $line . "\n" if $printing;
	}
	close $INALN;
	close $OUTALN;
	`mv "$alignment.seqids" "$alignment"`;
	foreach my $key ( keys %mapped_ids ) {
		$existing_ids->{$key} = $mapped_ids{$key};
	}
	return $counter;
}

sub read_gene_ids {
	my %args = @_;
	my $file = $args{file} || miss("file");
	my $IDS = ps_open( $file );
	my %id_to_taxon;
	my %marker_taxon_to_id;
	while ( my $line = <$IDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		$id_to_taxon{$uniqueid} = $taxon;
		push( @{ $marker_taxon_to_id{$marker}{$taxon} }, $uniqueid );
	}
	return ( \%id_to_taxon, \%marker_taxon_to_id );
}

sub update_ncbi_taxonomy {
	my %args = @_;
	my $repository = $args{repository} || miss("repository");
	my $base_marker_dir = $args{base_marker_dir} || miss("base_marker_dir");
	print "Downloading new NCBI taxonomy...\n";
	my $url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
#	print "Include path is @INC\n";
#	my $ff = File::Fetch->new( uri => $url );
	`mkdir -p "$repository/ncbi"`;
#	$ff->fetch( to => "$repository/ncbi" );
	`cd $repository/ncbi ; curl -LO $url`;
	system("cd $repository/ncbi/ ; tar xzf taxdump.tar.gz");
	unlink("$repository/ncbi/taxdump.tar.gz");
	unlink("$repository/ncbi/citations.dmp");
	unlink("$repository/ncbi/division.dmp");
	unlink("$repository/ncbi/gc.prt");
	unlink("$repository/ncbi/gencode.dmp");
	unlink("$repository/ncbi/readme.txt");
	system("cd $repository/; tar czf ncbi.tar.gz ncbi");
	`cd $base_marker_dir ; rm taxdmp.zip ; taxit new_database`;
}

sub make_marker_taxon_map {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $markerdir   = $args{marker_dir} || miss("marker_dir");

	my $AAIDS = ps_open( "$markerdir/".get_gene_id_file(dna=>0) );
	my $MARKERTAXONMAP = ps_open( ">$markerdir/marker_taxon_map.updated.txt" );

	while ( my $line = <$AAIDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		print $MARKERTAXONMAP "$uniqueid\t$taxon\n";
	}
	close $MARKERTAXONMAP;
}

sub make_ncbi_tree_from_update {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $markerdir   = $args{marker_dir} || miss("marker_dir");

	my $AAIDS = ps_open( "$markerdir/".get_gene_id_file(dna=>0) );
	my @taxonids;

	while ( my $line = <$AAIDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		push( @taxonids, $taxon ) if $taxon =~ /^\d+$/;
	}
	Phylosift::MarkerBuild::make_ncbi_subtree(out_file=>"$markerdir/ncbi_tree.updated.tre", taxon_ids=>\@taxonids);
}

sub get_marker_name_base {
	my %args = @_;
	my $name = $args{marker} || miss("marker");
	my $dna = $args{dna} || 0;
	my $updated = $args{updated} || 0;
	my $pruned = $args{pruned} || 0;
	my $sub_marker = $args{sub_marker};
	$name .= ".codon"   if $dna;
	$name .= ".updated" if $updated;

	# rna markers don't get pruned
	$name .= ".pruned" if $pruned && Phylosift::Utilities::is_protein_marker( marker => $args{marker} );
	$name .= ".sub".$sub_marker if defined($sub_marker);
	return $name;
}

sub get_marker_package {
	my %args = @_;
	return get_marker_name_base(%args);
}

sub get_marker_geneids {
	my %args = @_;
	my $name = get_marker_name_base(%args);
	$name .= ".geneids";
	return $name;
}

sub get_marker_ncbi_subtree {
	my %args = @_;
	my $name = get_marker_name_base(%args);
	$name .= ".ncbi_subtree";
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

sub build_marker_trees_fasttree {
	my %args        = @_;
	my $marker_dir  = $args{marker_directory} || miss("marker_directory");
	my $pruned      = $args{pruned} || 0;
	my $codon_fasta = get_fasta_filename( marker => '$1', dna => 1, updated => 1, pruned => 0, sub_marker => '$2' );
	my $aa_fasta    = get_fasta_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_tre   = get_fasttree_tre_filename( marker => '$1', dna => 1, updated => 1, pruned => 0, sub_marker => '$2' );
	my $aa_tre      = get_fasttree_tre_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );
	my $codon_log   = get_fasttree_log_filename( marker => '$1', dna => 1, updated => 1, pruned => 0, sub_marker => '$2' );
	my $aa_log      = get_fasttree_log_filename( marker => '$1', dna => 0, updated => 1, pruned => $pruned );	
	my $TREESCRIPT = ps_open( ">/tmp/ps_tree.sh" );
	print $TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
export OMP_NUM_THREADS=3
/home/koadman/bin/FastTree -log $aa_log $aa_fasta > $aa_tre
};
	`chmod 755 /tmp/ps_tree.sh`;
	close $TREESCRIPT;
	$TREESCRIPT = ps_open( ">/tmp/ps_tree_codon.sh" );
	print $TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
/home/koadman/bin/FastTree -nt -gtr -log $codon_log  $codon_fasta > $codon_tre
};
	`chmod 755 /tmp/ps_tree_codon.sh`;
	close $TREESCRIPT;
	$TREESCRIPT = ps_open( ">/tmp/ps_tree_rna.sh" );
	print $TREESCRIPT qq{#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
#\$ -pe threaded 3
export OMP_NUM_THREADS=3
/home/koadman/bin/FastTree -nt -gtr -log $aa_log $aa_fasta > $aa_tre
};
	`chmod 755 /tmp/ps_tree_rna.sh`;
	chdir($marker_dir);
	my @markerlist = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	my @jobids;

	foreach my $marker (@markerlist) {
		my $marker_fasta = get_fasta_filename( marker => $marker, dna => 0, updated => 1, pruned => $pruned );
		next unless ( -e $marker_fasta );
		my $tree_script = "/tmp/ps_tree.sh";
		$tree_script = "/tmp/ps_tree_rna.sh" unless Phylosift::Utilities::is_protein_marker( marker => $marker );

		# run fasttree on them
		my $qsub_args;
		$qsub_args = "-pe threaded 3" if $marker eq "concat";
		qsub_job(script=>$tree_script, job_ids=>\@jobids, qsub_args => $qsub_args, script_args=>[$marker] );
		next unless Phylosift::Utilities::is_protein_marker( marker => $marker );

		# only do codon markers once aa marker trees have already been pruned
		next unless $pruned;

		# run fasttree on codons
		for(my $group_id=1; ; $group_id++){
			my $subalignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			last unless -e $subalignment;
			my $mtime = ( stat($subalignment) )[9];
			last unless $mtime > 1332510000;
			my @sa = ($marker, $group_id);
			qsub_job(script=>"/tmp/ps_tree_codon.sh", job_ids=>\@jobids, script_args=>\@sa );
		}
			
	}
	wait_for_jobs( job_ids => \@jobids );
	`rm ps_tree.sh.* ps_tree_codon.* ps_tree_rna.*`;
}

sub create_taxon_id_table {
	my %args = @_;
	my $alignment = $args{alignment} || miss("alignment");
	my $output = $args{output} || miss("output");
	my $alltaxa = $args{alltaxa} || miss("alltaxa");
	my $OUTPUT = ps_open( ">$output" );
	my $INALN = ps_open( $alignment );
	while ( my $line = <$INALN> ) {
		if ( $line =~ /^>(.+)/ ) {
			my $header = $1;
			if ( $header =~ /\.(\d+?)\.fasta/ ) {
				print $OUTPUT "$header\t$1\n";
				$alltaxa->{$1}=1;
			}
		}
	}
	close $OUTPUT;
}

sub fix_names_in_alignment {
	my %args = @_;
	my $alignment = $args{alignment} || miss("alignment");
	my %markertaxa;

	# naive means to remove most trivial duplicates. TODO: allow divergent paralogs to remain.
	my $printing = 1;
	my $INALN = ps_open( $alignment );
	my $OUTALN = ps_open( ">$alignment.fixed" );
	while ( my $line = <$INALN> ) {
		if ( $line =~ /^>(.+)/ ) {
			my $header = $1;
			if ( $header =~ /\.(\d+?)\.fasta/ ) {
				$line           = ">$1\n";
				$printing       = defined( $markertaxa{$1} ) ? 0 : 1;
				$markertaxa{$1} = 1;
			}
		}
		print $OUTALN $line if $printing;
	}
	close $INALN;
	close $OUTALN;
	`mv "$alignment.fixed" "$alignment"`;
}

sub filter_fasta {
	my %args         = @_;
	my $input_fasta  = $args{input_fasta} || miss("input_fasta");
	my $output_fasta = $args{output_fasta} || miss("output_fasta");
	my $keep_taxa    = $args{keep_taxa} || miss("keep_taxa");

	# create a pruned fasta
	my $FASTA = ps_open( $input_fasta );
	my $PRUNEDFASTA = ps_open( ">" . $output_fasta );
	my $printing = 0;
	while ( my $line = <$FASTA> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			$printing = defined( $keep_taxa->{$1} ) ? 1 : 0;
		}
		print $PRUNEDFASTA "$line\n" if $printing;
	}
}

sub prune_marker {
	my %args         = @_;
	my $tre          = $args{tre} || miss("tre");
	my $distance     = $args{distance} || miss("distance");
	my $fasta        = $args{fasta} || miss("fasta");
	my $pruned_fasta = $args{pruned_fasta} || miss("pruned_fasta");

	# no point in pruning something with fewer than three taxa
	my $seq_count = `grep -c ">" "$fasta"`;
	chomp $seq_count;
	if($seq_count < 3){
		`cp "$fasta" "$pruned_fasta"`;
		return;
	}

	#	my $prune_cmd = "$Phylosift::Settings::pda -k 20000 -g -minlen $args{distance} $args{tre} $args{tre}.pruning.log";
	my $prune_cmd = "pda -k 20000 -g -minlen $distance $tre $tre.pruning.log";
	system("$prune_cmd");

	# read the list of taxa to keep
	my $PRUNE = ps_open( "$tre.pruning.log" );
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
	filter_fasta( input_fasta => $fasta, output_fasta => $pruned_fasta, keep_taxa => \%keep_taxa );
	`rm "$tre.pruning.log"`;
}

sub pd_prune_markers {
	my %args = @_;
	my $marker_dir = $args{marker_directory} || miss("marker_directory");
	chdir($marker_dir);

	# prune distance is different for AA and DNA
	# these distances are in substitutions per site
#	my %PRUNE_DISTANCE = ( 0 => 0.01, 1 => 0.003 );
	my %PRUNE_DISTANCE = ( 0 => 0.001, 1 => 0.003 );	# try leaving all but the totally identical taxa in to make jplace comparison possible
	my $REPS_DISTANCE  = 0.05;                                  # reps can be further diverged since we care only about similarity search and not read placement
	my @markerlist     = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	my $dna = 0;
#	for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {                   # zero for aa, one for dna
		foreach my $marker (@markerlist) {
			next unless Phylosift::Utilities::is_protein_marker( marker => $marker );
			my $tre = get_fasttree_tre_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
			my $fasta  = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );			
			my $pruned = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => 1 );
			next unless -e $fasta && -e $tre;	# perhaps not all markers have made it this far...
			
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
#	}
}

sub qsub_job {
	my %args = @_;
	my $script = $args{script} || miss("script");
	my $saref = $args{script_args};
	my $qsub_args = $args{qsub_args} || "";
	my $jobsref = $args{job_ids};
	
#	my $qsub_cmd = "qsub -q all.q -q eisen.q $qsub_args $script ";
	my $qsub_cmd = "qsub $qsub_args $script ";
	$qsub_cmd .= join(" ", @$saref ) if defined($saref);
	my $job      = `$qsub_cmd`;
	$job =~ /Your job (\d+) /;
	push( @$jobsref, $1 );
}

sub filter_marker_gene_ids{
	my %args = @_;
	my $marker = $args{marker} || miss("marker");
	my $updated = $args{updated} || 0;
	my $dna = $args{dna} || 0;
	my $pruned = $args{pruned} || 0;
	my $sub_marker = $args{sub_marker};
	
	my @taxon_ids;
	
	my $ids = get_marker_geneids( marker => $marker, dna => $dna, updated => $updated, sub_marker=>$sub_marker );
	print STDERR "ids are $ids\n";
	my $MARKER_IDS = ps_open(">$ids");

	# read the IDs present in the alignment so only those are included
	my $alignment = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned, sub_marker => $sub_marker );
	my $ALN = ps_open($alignment);
	my %aln_ids;
	while(my $line = <$ALN>){
		chomp $line;
		$aln_ids{$1}=1 if $line =~ />(.+)/;
	}
	my $merged = Phylosift::Summarize::read_merged_nodes();
	
	

	my $IDS = ps_open(get_gene_id_file(dna => $dna));
	while(my $line = <$IDS>){
		next unless $line =~ /^$marker\t/;
		chomp $line;
		my ( $m, $taxon, $uniqueid ) = split( /\t/, $line );
		next unless defined($aln_ids{$uniqueid});
		$taxon = $merged->{$taxon} if defined($merged->{$taxon});
		push(@taxon_ids, $taxon);
		print $MARKER_IDS "$m\t$taxon\t$uniqueid\n";
	}
	close($MARKER_IDS);	
	return \@taxon_ids;
}

sub reconcile_with_ncbi {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $marker_dir  = $args{marker_directory} || miss("marker_directory");
	my $pruned      = $args{pruned} || miss("pruned");
	
	my $p_orig = $pruned;
	chdir $marker_dir;
	my $aa_script   = "/tmp/ps_reconcile.sh";
	my $codon_script = "/tmp/ps_reconcile_codon.sh";
	for(my $dna=0; $dna<2; $dna++){
		my $script = $dna == 0 ? $aa_script : $codon_script;
		my $sub_marker;
		$sub_marker = '$2' if $dna;
		$pruned = 0 if $dna;
		my $aa_fasta    = get_fasta_filename( marker => '$1', dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
		my $aa_tre      = get_fasttree_tre_filename( marker => '$1', dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
		my $aa_log      = get_fasttree_log_filename( marker => '$1', dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
		my $aa_package  = get_marker_package( marker => '$1', dna => $dna, updated => 1, sub_marker=>$sub_marker );
		my $aa_taxonmap = get_taxonmap_filename( marker => '$1', dna => $dna, updated => 1, sub_marker=>$sub_marker );
		my $aa_ids      = get_marker_geneids( marker => '$1', dna => $dna, updated => 1, sub_marker=>$sub_marker );
		my $aa_tmpread  = get_marker_package( marker => '$1', dna => $dna, updated => 1, sub_marker=>$sub_marker ).".tmpread";
		my $aa_ncbi_tre = get_marker_ncbi_subtree( marker => '$1', dna => $dna, updated => 1, sub_marker=>$sub_marker );
		my $RECONCILESCRIPT = ps_open( ">$script" );
		print $RECONCILESCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

rm -rf $aa_package
taxit create -a "Aaron Darling" -d "simple package for reconciliation only" -l \$1 -f $aa_fasta -t $aa_tre -s $aa_log -P $aa_package
pplacer -c $aa_package -p $aa_tmpread.fasta
# readconciler uses a pplacer tree from a .jplace file to parse out the branch numbers
mangler.pl < $aa_tmpread.jplace > $aa_tmpread.jplace.mangled
/home/koadman/development/PhyloSift/tools/make_ncbi_subtree.pl $aa_ncbi_tre $marker_dir \$1 1 $pruned $dna $sub_marker
readconciler $aa_ncbi_tre $aa_tmpread.jplace.mangled $aa_ids $aa_taxonmap
# rm $aa_tmpread.jplace $aa_tmpread.jplace.mangled $aa_tmpread.fasta
EOF
		close($RECONCILESCRIPT);
		`chmod 755 $script`;
	}
	$pruned = $p_orig;

	my @markerlist = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	my @jobids;

	foreach my $marker (@markerlist) {
		my $marker_fasta = get_fasta_filename( marker => $marker, updated => 1, pruned=>$pruned);
		print STDERR "Couldnt find $marker_fasta\n" unless -e $marker_fasta;
		next unless -e $marker_fasta;
		next if $marker =~ /^PMPROK/ || $marker =~ /^DNGNGWU/;	# skip these since they all go in the concat
		my $taxa = filter_marker_gene_ids(marker=>$marker, updated=>1, pruned=>$pruned, dna=>0);
		print STDERR "working on marker $marker\n";
		# create some read files for pplacer to place so we can get its jplace
		Phylosift::MarkerBuild::create_temp_read_fasta( file => get_marker_package( marker=>$marker, updated=>1), aln_file => $marker_fasta);


		# run reconciliation on them
		my @marray = ($marker);
		my $qsub_args = $marker eq "concat" ? "-l mem_free=20G" : "";
		qsub_job(script=>$aa_script, qsub_args => $qsub_args, job_ids=>\@jobids, script_args=>\@marray );

		for(my $group_id=1; ; $group_id++){
			my $subalignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			last unless -e $subalignment;
			my $mtime = ( stat($subalignment) )[9];
			last unless $mtime > 1332510000;

			Phylosift::MarkerBuild::create_temp_read_fasta( file => get_marker_package( marker=>$marker, updated=>1, dna=>1, sub_marker=>$group_id), aln_file => get_fasta_filename( marker => $marker, updated => 1, dna => 1, sub_marker=>$group_id) );
			my $taxa = filter_marker_gene_ids(marker=>$marker, updated=>1, dna=>1, pruned=>0, sub_marker=>$group_id);
			my $ncbi_tre = get_marker_ncbi_subtree( marker => $marker, dna => 1, updated => 1, pruned=>0, sub_marker=>$group_id );
			my $success = -f $ncbi_tre && -s $ncbi_tre > 0;
			debug "Submitting reconciliation job for $marker codon sub $group_id\n";
			my @marray = ($marker,$group_id);
			qsub_job(script=>$codon_script, job_ids=>\@jobids, script_args=>\@marray );
		}
	}
	wait_for_jobs( job_ids => \@jobids );
	`rm ps_reconcile.sh.o* ps_reconcile.sh.e* ps_reconcile_codon.sh.o* ps_reconcile_codon.sh.e* *.tmpread.fasta *.tmpread.jplace *.tmpread.jplace.mangled`;
}

sub write_marker_build_script {
	my %args = @_;
	my $bs = $args{batch_script};
	my $BUILDSCRIPT = ps_open( ">$bs" );
	my $ps = $FindBin::Bin;
	print $BUILDSCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

export PATH="\$PATH:$ps"
$ps/phylosift build_marker -f --alignment=\$1 --update-only --taxonmap=\$2 --reps_pd=\$3 \$4

EOF

}

sub launch_marker_build {
	my %args = @_;
	my $dna = $args{dna};
	my $marker = $args{marker} || miss("marker");
	my $bs = $args{batch_script} || miss("batch_script");
		
	my @jobids;

	my $marker_fasta = get_fasta_filename( marker => $marker, updated => 1, dna=>$dna);
	print STDERR "Couldnt find $marker_fasta\n" unless -e $marker_fasta;
	next unless -e $marker_fasta;
	next if $marker =~ /^PMPROK/ || $marker =~ /^DNGNGWU/;	# skip these since they all go in the concat
	
	my $qsub_args = $marker eq "concat" ? "-l mem_free=20G" : "";
	my @marray = ($marker_fasta,$marker_fasta.".taxon_ids","0.01");
	my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
	push(@marray, "--unaligned=$clean_reps") if $dna==0;
	qsub_job(script=>$bs, qsub_args => $qsub_args, job_ids=>\@jobids, script_args=>\@marray );
	return $jobids[0];
}

sub launch_marker_builds {
	my %args = @_;
	my $self        = $args{self} || miss("self");
	my $marker_dir  = $args{marker_dir} || miss("marker_dir");
		
	chdir("$marker_dir");
	my @markerlist = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	my @jobids;

	my $bs = "/tmp/ps_build_marker.sh";
	my $BUILDSCRIPT = ps_open( ">$bs" );
	my $ps = $FindBin::Bin;
	print $BUILDSCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

$ps/phylosift build_marker -f --alignment=\$1 --update-only --taxonmap=\$2 --reps_pd=\$3 \$4

EOF

	for(my $dna=0; $dna<2; $dna++){
		foreach my $marker (@markerlist) {
			my $marker_fasta = get_fasta_filename( marker => $marker, updated => 1, dna=>$dna);
			print STDERR "Couldnt find $marker_fasta\n" unless -e $marker_fasta;
			next unless -e $marker_fasta;
			next if $marker =~ /^PMPROK/ || $marker =~ /^DNGNGWU/;	# skip these since they all go in the concat
			
			my $qsub_args = $marker eq "concat" ? "-l mem_free=20G" : "";
			my @marray = ($marker_fasta,$marker_fasta.".taxon_ids","0.01");
			my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
			push(@marray, "--unaligned=$clean_reps") if $dna==0;
			qsub_job(script=>$bs, qsub_args => $qsub_args, job_ids=>\@jobids, script_args=>\@marray );
		}
	}
	wait_for_jobs( job_ids => \@jobids );
}

=head2 update_rna

RNA markers contain sequences for many taxa that aren't available in genome databases.
We want the updated markers to retain these original sequences.

=cut

sub update_rna {
	my %args        = @_;
	my $self        = $args{self} || miss("self");
	my $marker_dir  = $args{marker_dir} || miss("marker_dir");
	my @marker_list = Phylosift::Utilities::gather_markers();
	foreach my $marker (@marker_list) {
		next if Phylosift::Utilities::is_protein_marker( marker => $marker );
		my $base_alignment = Phylosift::Utilities::get_marker_aln_file( self=> $self, marker=>$marker );

		#		my $base_alignment = Phylosift::Utilities::get_marker_aln_file(self=>$self, marker=>$marker);
		my $updated_alignment = get_fasta_filename( marker=>$marker, updated=>1 );

		# need to add to gene_ids.aa.txt
		# use greengenes id for taxon ID
		my $aaid_file = "$marker_dir/".get_gene_id_file(dna=>0);
		my $max_id    = -1;

		#
		# get the largest gene ID
		my $AAID = ps_open( $aaid_file );
		while ( my $line = <$AAID> ) {
			chomp $line;
			my ( $marker, $tid, $gid ) = split( /\t/, $line );
			$max_id = $gid if $gid > $max_id;
		}
		close $AAID;

		#
		# now start assigning IDs to base marker sequences and adding them to the update
		$max_id++;
		$AAID = ps_open( ">>" . $aaid_file );
		my $INALN = ps_open( $base_alignment );
		my $OUTALN = ps_open( ">>" . $updated_alignment );
		while ( my $line = <$INALN> ) {
			chomp $line;
			if ( $line =~ /^>(.+)/ ) {
				my $countstring = sprintf( "%010u", $max_id );
				my $ggid = $1;
				$ggid =~ s/\/\d+-\d+//g;    # remove trailing bioperl rubbish
				print $AAID "$marker\tgg_$ggid\t$countstring\n";
				$line = ">$countstring";
				$max_id++;
			}
			print $OUTALN "$line\n";
		}
		close $AAID;
	}
}

=head2 make_codon_submarkers

Find groups of closely related taxa that would be better analyzed with DNA sequence rather than protein

=cut

sub make_codon_submarkers {
	my %args = @_;
	my $marker_dir = $args{marker_dir} || miss("marker_dir");
	chdir($marker_dir);

	# this is the maximum distance in amino acid substitutions per site that are allowed
	# on any branch of the tree relating members of a group
	# TODO: tune this value
	my $max_aa_branch_distance = 0.2;
	my ($idtref, $mtiref) = read_gene_ids( file => "$marker_dir/".get_gene_id_file(dna=>0) );
	my %id_to_taxon = %$idtref;
	my %marker_taxon_to_id = %$mtiref;

	($idtref, $mtiref) = read_gene_ids( file => "$marker_dir/".get_gene_id_file(dna=>1) );
	my %codon_id_to_taxon = %$idtref;
	my %codon_marker_taxon_to_id = %$mtiref;

	my @subalignments;    # list of subalignments created that will later need marker packages made
	
	# the following file will provide a mapping of AA gene ID to submarker
	my $SUBTABLE = ps_open( ">submarkers.txt" );
	my @marker_list = Phylosift::Utilities::gather_markers();
	unshift( @marker_list, "concat" );
	foreach my $marker (@marker_list) {
		next unless Phylosift::Utilities::is_protein_marker( marker => $marker );
		my $aa_tree = get_fasttree_tre_filename( marker => $marker, dna => 0, updated => 1, pruned => 0 );
		my $codon_alignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0 );
		# clean out any stale version of this marker
		my $sub_pack = get_marker_package( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker=>'*' );
		debug "removing $sub_pack\n";
		`rm -rf $sub_pack`;
		next unless -e $codon_alignment;
		my $alnio = Bio::AlignIO->new(-file => $codon_alignment );
		my $aln = $alnio->next_aln();

		# get groups of taxa that are close
		my $GROUP_TABLE = ps_open("segment_tree $aa_tree $max_aa_branch_distance |");
		my %gene_groups;
		my $group_id = 1;
		while( my $group_line = <$GROUP_TABLE>) {			
			chomp($group_line);
			my @gene_ids = split( /\t/, $group_line );

			# create a new subalignment file
			my $subalignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			my $SUBALN = ps_open( ">$subalignment" );

			# map the gene ID from aa tree into the corresponding ID in the codon data
			foreach my $gene (@gene_ids) {
				my $taxon     = $id_to_taxon{$gene};
				if(!defined($codon_marker_taxon_to_id{$marker}{$taxon})){
					croak("Error: codon table missing taxon $taxon in marker $marker");
				}
				my @codon_ids = @{ $codon_marker_taxon_to_id{$marker}{$taxon} };

				# write each sequence into the subalignment
				foreach my $id (@codon_ids) {
					foreach my $seq ( $aln->each_seq_with_id($id) ) {
						print $SUBALN ">$id\n";
						print $SUBALN $seq->seq() . "\n";
					}
					my $codon_taxon     = $codon_id_to_taxon{$id};
					if($codon_taxon ne $taxon){
						croak("Error: inconsistent taxa in AA and Codon data. AA: $taxon, Codon: $codon_taxon, aa gene: $gene, codon gene: $id, marker: $marker");
					}
				}
				

				# add the mapping from gene ID to submarker to the table
				print $SUBTABLE "$gene\t$marker\t$taxon\t$group_id\n";
			}
			close $SUBALN;
			$group_id++;
		}
	}
	close $SUBTABLE;
}

=head2 make_constrained_tree

Makes a tree for a marker that respects topological constraints given by another marker's tree

=cut

sub make_constrained_tree {
	my %args              = @_;
	my $constraint_marker = $args{constraint_marker} || miss("constraint_marker");
	my $constraint_tree   = $args{constraint_tree} || miss("constraint_tree");
	my $target_marker     = $args{target_marker} || miss("target_marker");
	my $marker_dir        = $args{marker_dir} || miss("marker_dir");
	my $constraint_pruned = $args{constraint_pruned} || 0;
	my $target_pruned     = $args{target_pruned} || 0;
	my $copy_tree         = $args{copy_tree} || 0;
	my ($idtref,$mtiref) = read_gene_ids(file=>"$marker_dir/".get_gene_id_file(dna=>0));
	my %id_to_taxon = %$idtref;
	my %marker_taxon_to_id = %$mtiref;

	croak("$constraint_tree does not exist") unless -e $constraint_tree;
	# first convert the constraint tree leaf names to match the target
	my $tree = Bio::Phylo::IO->parse(
									  '-file'   => $constraint_tree,
									  '-format' => 'newick',
	)->first;
	print STDERR "Tree has ".scalar(@{ $tree->get_entities })." nodes\n";
	
	# collect all the sequence IDs present in the target alignment
	my %target_ids;
	my $target_alignment         = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned );
	my $target_alignment_amended = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ) . ".amended";
	my $column_count             = 0;
	my $TALN = ps_open($target_alignment);
	my $dat = "";
	while(my $line = <$TALN>){
		chomp $line;
		if($line =~ /^>(.+)/){
			$target_ids{$1}=1;
			$column_count = length($dat) if($column_count == 0);
		}elsif($column_count==0){
			$dat .= $line;
		}
	}
	
	my @missing_taxa;    # gene IDs for taxa missing in the target
	foreach my $node ( @{ $tree->get_entities } ) {
		# skip this one if it is not a leaf
		next if ( scalar( @{ $node->get_children() } ) > 1 );
		my $name = $node->get_name;
		croak "Unable to find taxon for gene $name in marker $constraint_marker" unless defined $id_to_taxon{$name};
		my $taxon_id = $id_to_taxon{$name};

		# now find the taxon in the target data, if it exists
		my $rename;
		if ( defined( $marker_taxon_to_id{$target_marker}{$taxon_id} ) ) {
			my @tids = @{ $marker_taxon_to_id{$target_marker}{$taxon_id} };
			$rename = $tids[0];
			# if it does exist, and the target is pruned, check whether it was pruned out
			push( @missing_taxa, $rename ) unless defined($target_ids{$rename});
		} else {
			push( @missing_taxa, $name );
		}
		$node->set_name($rename);
	}
	my $renamed_tree = "$constraint_tree.renamed";
	my $RENAMETREE = ps_open( ">$renamed_tree" );
	print $RENAMETREE unparse( '-phylo' => $tree, '-format' => 'newick' ) . "\n";
	close $RENAMETREE;

	
	#
	# create an amended alignment with any taxa missing from the target data
	`cp "$target_alignment" "$target_alignment_amended"`;
	my $AMENDALN = ps_open( ">>$target_alignment_amended" );
	foreach my $missing (@missing_taxa) {
		print $AMENDALN ">$missing\n" . ( "-" x $column_count ) . "\n";
	}
	close $AMENDALN;

	unless($copy_tree){
		#
		# enumerate the splits in the constraint tree
		my $constraint_splits = "$renamed_tree.splits";
		my $constraint_cl     = "print_splits $renamed_tree $constraint_splits";
		system($constraint_cl);
		#
		# now infer a target tree that respects protein tree constraints
		my $target_constrained_tree = get_fasttree_tre_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ) . ".constrained";
		my $target_constrained_log = get_fasttree_log_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ) . ".constrained";
		
		my $data_type_args = "";
		$data_type_args = "-nt -gtr" unless(Phylosift::Utilities::is_protein_marker(marker=>$target_marker));
		
		my $TREESCRIPT = ps_open(">/tmp/constrained_tre.sh");
		print $TREESCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash
#\$ -pe threaded 3
export OMP_NUM_THREADS=3
FastTree-nonuni $data_type_args -constraints $constraint_splits -constraintWeight 100 -log $target_constrained_log  $target_alignment_amended > $target_constrained_tree
	
EOF
		my @job_ids = ();
		qsub_job(script=>"/tmp/constrained_tre.sh", job_ids=>\@job_ids, qsub_args=>"-l mem_free=12G");
		wait_for_jobs(job_ids=>\@job_ids);

		# finally, make a pplacer package with the constrained tree
		my $taxit_cl =
	"taxit create -a \"Aaron Darling\" -d \"topology-constrained marker $target_marker\" -l $target_marker -f $target_alignment_amended -t $target_constrained_tree -s $target_constrained_log -P $target_marker.constrained";
		system($taxit_cl);
		unlink($constraint_splits);
	}else{

		# make a pplacer package with the renamed tree
		my $taxit_cl =
	"taxit create -a \"Aaron Darling\" -d \"topology-constrained marker $target_marker\" -l $target_marker -f $target_alignment_amended -t $renamed_tree -P $target_marker.constrained";
		system($taxit_cl);		
	}
}

=head2 join_trees

Create protein & rna trees with compatible topologies

=cut

sub join_trees {
	my %args = @_;
	my $marker_dir = $args{marker_dir} || miss("marker_dir");

	# first apply protein constraints to 16s
	# then apply 16s constraints to protein
#	my $constraint_tree = get_fasttree_tre_filename( marker => "concat", dna => 0, updated => 1, pruned => 1 );
#	make_constrained_tree( constraint_marker => "concat", constraint_tree => $constraint_tree,   target_marker => "16s_reps_bac", marker_dir => $marker_dir, constraint_pruned => 1, target_pruned => 0 );
#	my $tctree = get_fasttree_tre_filename( marker => "16s_reps_bac", dna => 0, updated => 1, pruned => 0 ) . ".constrained";
	my $tctree = get_fasttree_tre_filename( marker => "16s_reps_bac", dna => 0, updated => 1, pruned => 0 );
	print "Using $tctree to provide constraints for concat\n";
	make_constrained_tree( constraint_marker => "16s_reps_bac", constraint_tree => $tctree, target_marker => "concat",       marker_dir => $marker_dir, constraint_pruned => 0, target_pruned => 1, copy_tree=>1 );

	# with a little luck we can now??? do what?
}

sub get_gene_id_file {
	my %args = @_;
	my $dna = $args{dna} || 0;
	return "gene_ids.codon.txt" if $dna;
	return "gene_ids.aa.txt";
}

sub make_pplacer_package {
	my %args = @_;
	my $dna = $args{dna} || 0;
	my $seq_ids = $args{seq_ids}; # optional
	my $marker = $args{marker} || miss("marker");
	my $pruned = $args{pruned} || 0;
	my $sub_marker = $args{sub_marker};

	print STDERR "packaging $marker sub $sub_marker dna $dna\n";
	# create a taxon id list for this marker
	my $gene_id_file = get_gene_id_file( dna => $dna );
	my $AAIDS = ps_open(  $gene_id_file );
	my $TAXIDS = ps_open( ">tax_ids.txt" );
	my $SEQINFO = ps_open( ">seq_info.csv" );
	while ( my $line = <$AAIDS> ) {
		chomp $line;
		my @dat = split( /\t/, $line );
		next unless $dat[0] eq $marker;
		print STDERR "found concat\n";
		if(defined($seq_ids)){
			next unless defined($seq_ids->{$dat[2]});
		}
		print $TAXIDS $dat[1] . "\n";
		print $SEQINFO "$dat[2],$dat[0],$dat[1],\n";
	}
	close $TAXIDS;
	close $SEQINFO;
	my $taxtable_cl = "taxit taxtable -d taxonomy.db -t tax_ids.txt -o taxa.csv";
	print STDERR "$taxtable_cl\n";
	system($taxtable_cl);

	# gather filenames to stuff into marker package
	my $fasta = get_fasta_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
	my $tre = get_fasttree_tre_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
	my $log = get_fasttree_log_filename( marker => $marker, dna => $dna, updated => 1, pruned => $pruned, sub_marker=>$sub_marker );
	my $pack = get_marker_package( marker => $marker, dna => $dna, updated => 1, sub_marker=>$sub_marker );
	my $taxit_cl =
"taxit create -l $marker -P $pack --taxonomy taxa.csv --seq-info seq_info.csv --tree-stats $log --tree-file $tre --aln-fasta $fasta ";
	print STDERR "$taxit_cl\n";
	system($taxit_cl);
}

sub make_pplacer_packages_with_taxonomy {
	my %args       = @_;
	my $marker_dir = $args{marker_dir} || miss("marker_dir");
	my $pruned     = $args{pruned} || 1;

	# only need to run this if it hasnt been done already
	# -- or do we need to re-do it for every taxonomy update?
	#	my $newdb_cl = "taxit new_database -d taxonomy.db";
	#	system($newdb_cl);
	my @markerlist = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	foreach my $marker (@markerlist) {
		# make the aa package
		make_pplacer_package(marker=>$marker, pruned=>$pruned);
		# now make any codon subpackages
		for(my $group_id=1; ; $group_id++){
			my $subalignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			last unless -e $subalignment;
			
			# extract the sequence IDs present in this submarker
			my %seq_ids;
			my $SUBALN = ps_open($subalignment);
			while( my $line = <$SUBALN>){
				chomp $line;
				if($line =~ />(.+)/){
					$seq_ids{$1}=1;
				}
			}
			make_pplacer_package(dna=>1, marker=>$marker, seq_ids=>\%seq_ids, pruned=>0, sub_marker => $group_id);
		}
	}
}

sub package_markers {
	my %args = @_;
	my $marker_dir = $args{marker_directory}|| miss("marker_directory");
	my $base_marker_dir = $args{base_marker_directory}|| miss("base_marker_directory");
	my @markerlist = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );
	
	# first collect the base markers
	foreach my $marker (@markerlist) {
		# copy the base marker into our directory 
		`cp -r $base_marker_dir/$marker/ $marker_dir`;
#		if($marker =~ /PMPROK/){
#			`cp $base_marker_dir/$marker.faa $marker_dir`;
#			`cp $base_marker_dir/$marker.ali $marker_dir`;
#			`cp $base_marker_dir/$marker.final.tre $marker_dir`;
#			`cp $base_marker_dir/$marker.ncbimap $marker_dir`;
#			`cp $base_marker_dir/$marker.stk $marker_dir`;
#			`cp $base_marker_dir/$marker.hmm $marker_dir`;
#			`cp $base_marker_dir/$marker.trimfinal $marker_dir`;
#			`cp $base_marker_dir/$marker.trimfinal.fasta $marker_dir`;
#			`cp $base_marker_dir/$marker.stock.hmm $marker_dir`;
#			`cp $base_marker_dir/$marker.seed.stock $marker_dir`;
#		}
	}

	# then move various files into place
	foreach my $marker (@markerlist) {
		for ( my $dna = 0 ; $dna < 2 ; $dna++ ) {    # zero for aa, one for dna
			my $pack = get_marker_package( marker => $marker, dna => $dna, updated => 1 );
			# move in reps
			my $pruned_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1, pruned => 1 );
			`mv "$pruned_reps" "$pack/$marker.reps"` if $dna==0 && -e $pruned_reps;
			# move in taxonmap
			my $taxonmap = get_taxonmap_filename( marker => $marker, dna => $dna, updated => 1, pruned => 0 );
			`mv "$taxonmap" "$pack/$marker.taxonmap"` if -e $taxonmap;
		}
		# now do codon submarkers
		for(my $group_id=1; ; $group_id++){
			my $subalignment = get_fasta_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			my $pack = get_marker_package( marker => $marker, dna => 1, updated => 1, sub_marker => $group_id );
			if(-d $pack && !-f "$pack/$subalignment"){
				print STDERR "Removing package $pack";
				`rm -rf "$pack"*`;
			}
			last unless -e $subalignment;
			# move in taxonmap
			my $taxonmap = get_taxonmap_filename( marker => $marker, dna => 1, updated => 1, pruned => 0, sub_marker => $group_id );
			`mv "$taxonmap" "$pack/$marker.taxonmap"` if -e $taxonmap;
		}

	}
	
	# clean up some stuff that really shouldn't be packaged
	system("rm -f *.fasttree.log");
	system("rm -f *.updated*fasta");
	system("rm -f *.reps.clean");
	system("rm -f *.updated.reps");
	system("rm -f *.updated*geneids");
	system("rm -f core.*");
	system("rm -f *.tmpread.jplace");
	system("rm -f *.ncbi_subtree");
	system("mv ncbi_tree.updated.tre ncbi.safekeeping");
	system("mv concat.updated.tre concat.safekeeping");
	system("rm -f *.updated*tre");
	system("mv concat.safekeeping concat.updated.tre");
	system("mv ncbi.safekeeping ncbi_tree.updated.tre");
	
	chdir( $marker_dir . "/../" );
	system("pwd");
	my @timerval = localtime();
	my $datestr  = Phylosift::Utilities::get_date_YYYYMMDD;
	system("tar czf markers_$datestr.tgz markers");
	`rm -f markers.tgz`;
	`ln -s markers_$datestr.tgz markers.tgz`;
}
1;
