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
use FindBin;
use Scalar::Util qw(looks_like_number);

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::UpdateDB - Functionality to download new genomes and update the marker database

=head1 VERSION

Version 0.01

=cut

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

	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/organelle.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/virus.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/phage.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaea.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/archaealvirus.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/bacteria.details.txt" );
	get_ebi_from_list( url => "http://www.ebi.ac.uk/genomes/eukaryota.details.txt" );
}

sub get_ebi_from_list {
	my %args = @_;
	my $list_url = $args{url} || miss("url");    # URL to EBI's table of genome characteristics
	`wget $list_url -O "ebi.list.txt"`;
	my $DETAILS = ps_open("ebi.list.txt");
	my $line    = <$DETAILS>;
	my %taxa_details;    # {$taxon_id}{$accession} = [$date, $version, $size]
	while ( $line = <$DETAILS> ) {

		# each line in the list file records a genome and some metadata
		chomp $line;
		my ( $acc, $version, $date, $taxid, $description ) = split( /\t/, $line );
		next if $description =~ /itochondri/;    # skip mitochondrion
		next if $description =~ /hloroplas/;     # skip chloroplast
		$taxa_details{$taxid}{$acc} = { dater => $date, version => $version };
	}
	eval {
		require LWP::Simple;
		LWP::Simple->import();
	};
	my $result = $@;
	die "Unable to load LWP::Simple" unless $result;
	foreach my $taxid ( keys(%taxa_details) ) {

		# add up the accession sizes
		my $total_size = 0;
		my $facc;
		my $date;

		foreach my $acc ( keys( %{ $taxa_details{$taxid} } ) ) {

			# if it's too big, don't download
			my $url      = "http://www.ebi.ac.uk/ena/data/view/$acc&display=txt&expanded=true";
			my $ua       = LWP::UserAgent->new( max_size => 200 );
			my $response = $ua->get($url);
			my $content  = $response->decoded_content();
			my $sizer    = $1 if $content =~ /(\d+) BP/;
			$total_size += $sizer;
			$facc = $acc unless defined $facc;
			print STDERR "taxid $taxid, acc $acc\n";
			$date = $taxa_details{$taxid}{$acc}->{dater} unless defined $date;
		}
		print STDERR "Total genome size of taxid $taxid is $total_size\n";

		next if $total_size > 200000000;    # 100Mbp limit

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
		foreach my $acc ( keys( %{ $taxa_details{$taxid} } ) ) {
			my $WGETTER = ps_open("| wget -i - -O $outfile.embl.tmp ");
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
	`rm "ebi.list.txt"`;
}

sub get_taxid_from_gbk {
	my %args = @_;
	my $file = $args{file} || miss("file");

	# first get the taxon ID
	my $GBK = ps_open($file);
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
	debug "CHANIGNG DIR : $directory\n";

	# First download all finished bacterial genomes
	# then for each genome, concatenate all replicons into a FastA file
	# Name the FastA file with the organism and its taxon ID.
	chdir($directory);
	my $ncbi_wget_cmd = "wget -m --continue --timeout=20 --accept=gbk ftp://ftp.ncbi.nih.gov/genomes/Bacteria";

	#		`$ncbi_wget_cmd`;
	my $FINDER = ps_open("find ftp.ncbi.nih.gov/genomes/Bacteria -type d |");
	while ( my $line = <$FINDER> ) {
		chomp $line;
		next unless $line =~ /\/Bacteria\/(.+)/;    #skipping if there is nothing after Bacteria
		my $orgname = $1;
		next if ( $orgname =~ /\// );               # could be a malformed genbank directory
		next unless length($1) > 1;
		my $seq_out;
		my $fasta_name;
		my $LSSER = ps_open("ls $line/*gbk |");

		while ( my $gbk = <$LSSER> ) {
			chomp $gbk;
			debug "Looking at $gbk\n";
			if ( !defined($fasta_name) ) {
				debug "constructing fasta_name\n";
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
			if ( -s "$gbk" ) {
				debug "Opening $gbk\n";
				my $seq_in = Phylosift::Utilities::open_SeqIO_object( file => "$gbk", format => 'GENBANK' );
				while ( my $inseq = $seq_in->next_seq ) {
					$seq_out->write_seq($inseq);
				}
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

	`$ncbi_wget_cmd`;
	my $FINDER = ps_open("find ftp.ncbi.nih.gov/genbank/wgs -name \"*.gz\" |");
	while ( my $line = <$FINDER> ) {
		next if $line =~ /\.mstr\.gbff\.gz/;
		chomp $line;
		`gunzip -f $line`;
		$line =~ s/\.gz//g;
		my $orgname = $1 if $line =~ /wgs\.(....)\./;
		my $taxid = get_taxid_from_gbk( file => $line );
		my $fasta_name = "$orgname.$taxid.fasta";
		if ( -e $fasta_name && ( stat($fasta_name) )[9] > ( stat($line) )[9] ) {
			print STDERR "Already have $fasta_name\n";
			last;
		}
		my $seq_out = Phylosift::Utilities::open_SeqIO_object( file => ">$fasta_name", format => "FASTA" );
		my $seq_in  = Phylosift::Utilities::open_SeqIO_object( file => "$line",        format => "genbank" );
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
	my $FINDER = ps_open("find . -name \"*.gbk\" |");

	while ( my $line = <$FINDER> ) {
		chomp $line;
		my $taxid = get_taxid_from_gbk( file => $line );

		# then unpack all the files and combine them
		my $fasta_out = basename( $line, ".gbk" ).".$taxid.fasta";
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
		`cat $catline >> "$fasta_out"` if ( $tarstatus == 0 );
		`rm $catline`;
	}
}

=head2 set_default_values

Sets default values for UpdateDB parameters

=cut

sub set_default_values {
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::MAX_SGE_JOBS,    value => 5000 );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::WORKDIR_normal,  value => "/state/partition1/koadman/phylosift/" );
	Phylosift::Settings::set_default( parameter => \$Phylosift::Settings::WORKDIR_fatnode, value => "/data/scratch/koadman/phylosift/" );

}

sub find_new_genomes {
	my %args        = @_;
	my $genome_dir  = $args{genome_directory} || miss("genome_directory");
	my $results_dir = $args{results_directory} || miss("results_directory");
	my $files       = $args{files} || miss("files");
	my $FINDER      = ps_open("find $genome_dir |");
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

sub qsub_updates {
	my %args            = @_;
	my $local_directory = $args{local_directory} || miss("local_directory");
	my $files           = $args{files} || miss("files");
	my $extended        = $args{extended} || 0;
	my $params          = "";

	#$params = " --updated_markers=0 --extended " if $extended;

	$params = " --isolate --besthit " if !$extended;

	# qsub a slew of phylosift jobs and have them report results to a local directory on this host for concatenation
	my $hostname = `hostname`;
	chomp $hostname;
	`mkdir -p "$local_directory"`;

	my %jobids;
	my $PHYLOSIFTSCRIPT = ps_open(">/tmp/pssge.sh");
	my $ps              = $FindBin::Bin;
	print $PHYLOSIFTSCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

WORKDIR=$Phylosift::Settings::WORKDIR_normal\$JOB_ID
mkdir -p \$WORKDIR
if [ \$? -gt 0 ]; then
   WORKDIR=$Phylosift::Settings::WORKDIR_fatnode\$JOB_ID
   mkdir -p \$WORKDIR
fi
cd \$WORKDIR

export PATH="\$PATH:$ps"
$ps/phylosift search $params \$1
$ps/phylosift align -f $params \$1
$ps/phylosift name $params \$1
rm -rf PS_temp/*/treeDir
rm -rf PS_temp/*/blastDir
rm -rf PS_temp/*/isolates.fasta
rm -rf PS_temp/*/alignDir/DNGNG*.newCandidate

scp -r PS_temp/* $hostname:$local_directory
#cp -r PS_temp/* $local_directory
rm -rf \$WORKDIR

EOF

	my $job_count = 0;

	foreach my $file ( @{$files} ) {
		$job_count++;
		qsub_job( script => "/tmp/pssge.sh", job_ids => \%jobids, script_args => [$file] );

		# check whether we've hit the limit for queued jobs, and rest if needed
		if ( $job_count == $Phylosift::Settings::MAX_SGE_JOBS ) {
			wait_for_jobs( job_ids => \%jobids, min_remaining => 20 );
			%jobids    = ();
			$job_count = 0;
		}
	}
	wait_for_jobs( job_ids => \%jobids );
}

sub wait_for_jobs {
	my %args          = @_;
	my $job_ids       = $args{job_ids} || miss("job_ids");
	my $min_remaining = $args{min_remaining} || 0;

	# wait for all (or most) jobs to complete
	while ( scalar( keys( %{$job_ids} ) ) > $min_remaining ) {
		foreach my $jobid ( keys( %{$job_ids} ) ) {
			my $output = `qstat -j $jobid 2>&1`;
			if ( !defined($output) || $output =~ /Following jobs do not exist/ ) {
				delete $job_ids->{$jobid};
			}
		}
		sleep(60);
	}
}

sub concat_marker_files {
	my %args            = @_;
	my $marker          = $args{marker};
	my $local_directory = $args{local_directory} || miss("local_directory");
	my $cfref           = $args{catfiles};
	my $cat_ch          = $args{cat_ch};

	my $catline = join( " ", @$cfref );
	$catline .= " ";
	my $fasta = get_fasta_filename( marker => $marker, updated => 1 );
	`cat $catline $cat_ch "$local_directory/$fasta" 2> /dev/null`;

	#print STDERR "\n\n\n\n\nnormal CATLINE\t $catline\n";
	$catline =~ s/updated\.(\d+)\.fasta /codon\.updated\.$1\.fasta /g;

	#print STDERR "\n\n\n\n\n\n\ncodon catline\t $catline\n";
	#exit;
	if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
		my $codon_fasta = get_fasta_filename( marker => $marker, updated => 1, dna => 1 );

		#print "CODON cat line\n $catline\n";
		`cat $catline $cat_ch "$local_directory/$codon_fasta" 2> /dev/null`;
	}
	unless ( $marker eq "concat" ) {
		$catline =~ s/\.codon\.updated\.fasta/\.unmasked/g;
		my $reps = get_reps_filename( marker => $marker, updated => 1 );
		`cat $catline $cat_ch "$local_directory/$reps" 2> /dev/null`;
	}
}

sub collate_markers {
	my %args            = @_;
	my $local_directory = $args{local_directory} || miss("local_directory");
	my $marker_dir      = $args{marker_dir} || miss("marker_dir");
	my $taxon_knockouts = $args{taxon_knockouts};

	# make the directory if it doesn't yet exist
	`mkdir -p $marker_dir`;

	# get list of markers
	chdir($marker_dir);
	my @markerlist = Phylosift::Utilities::gather_markers();
	print STDERR "Markers are ".join( " ", @markerlist )."\n";

	my %ko_list;
	if ( defined($taxon_knockouts) ) {
		my $KO = ps_open($taxon_knockouts);
		while ( my $line = <$KO> ) {
			chomp $line;
			$ko_list{$line} = 1;
		}
	}

	my $bs = "/tmp/ps_build_marker.sh";
	write_marker_build_script( batch_script => $bs, destination => $marker_dir );

	# get a list of genomes available in the results directory
	# this hopefully means we touch each inode over NFS only once
	# NFS is slooooow...
	print STDERR "Working with ".scalar(@markerlist)." markers\n";
	print STDERR "Listing all files in results dir\n";
	my @alldata = `find "$local_directory" -name "*.fasta"`;
	print STDERR "Found ".scalar(@alldata)." sequence files\n";
	unshift( @markerlist, "concat" );
	my %alltaxa;
	my %job_ids;

	foreach my $marker (@markerlist) {
		my $cat_ch = ">";    # first time through ensures that existing files get clobbered

		# find all alignments with this marker
		my @catfiles = ();
		foreach my $file (@alldata) {
			next unless $file =~ /(.+\.fasta)\/alignDir\/$marker.updated.\d+.fasta/;
			my $genome = $1;
			my $taxon = $1 if $genome =~ /\.(\d+)\.fasta/;
			next if ( defined($taxon) && defined( $ko_list{$taxon} ) );
			chomp($file);
			push( @catfiles, $file );

			# cat the files into the alignment in batches
			# this cuts down on process spawning and file I/O
			# can't do more than about 4k at once because argument lists get too long
			my $batch_size = 1000;
			if ( @catfiles > $batch_size ) {
				print STDERR "Found $batch_size files for marker $marker\n";
				concat_marker_files( marker => $marker, local_directory => $local_directory, catfiles => \@catfiles, cat_ch => $cat_ch );
				@catfiles = ();
				$cat_ch   = ">>";    # now add to existing files
			}
		}
		if ( @catfiles > 0 ) {
			print STDERR "Last cat for marker $marker\n";
			concat_marker_files( marker => $marker, local_directory => $local_directory, catfiles => \@catfiles, cat_ch => $cat_ch );
		}

		# now rename sequences with their taxon IDs
		my $fasta = get_fasta_filename( marker => $marker, updated => 1 );
		next unless -e "$local_directory/$fasta";
		if ( Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
			filter_short_and_unclassified_seqs_from_fasta( input_fasta => "$local_directory/$fasta", output_fasta => "$marker_dir/$fasta", min_pct => 25 );
		} else {
			`cp $local_directory/$fasta $marker_dir/$fasta`;
		}
		create_taxon_id_table( alignment => "$marker_dir/$fasta", output => "$marker_dir/$fasta.taxon_ids", alltaxa => \%alltaxa );

		if ( !Phylosift::Utilities::is_protein_marker( marker => $marker ) ) {
			update_rna( marker => $marker, marker_dir => $marker );
		}

		#		`cp "$local_directory/$fasta" "$marker_dir/$fasta"`;
		my $reps = get_reps_filename( marker => $marker, updated => 1 );
		my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
		if ( -e "$local_directory/$reps" ) {
			clean_representatives( infile => "$local_directory/$reps", outfile => "$local_directory/$clean_reps" );
			`cp "$local_directory/$clean_reps" "$marker_dir/$clean_reps"`;
		}
		debug "Launching marker build for $marker\n";
		my $job_id = launch_marker_build( marker => $marker, dna => 0, batch_script => $bs );
		$job_ids{$job_id} = 1 if defined $job_id;

		my $codon_fasta = get_fasta_filename( marker => $marker, updated => 1, dna => 1 );
		next unless -e "$local_directory/$codon_fasta";
		filter_short_and_unclassified_seqs_from_fasta(
													   input_fasta  => "$local_directory/$codon_fasta",
													   output_fasta => "$marker_dir/$codon_fasta",
													   min_pct      => 25
		);
		create_taxon_id_table( alignment => "$marker_dir/$codon_fasta", output => "$marker_dir/$codon_fasta.taxon_ids", alltaxa => \%alltaxa );

		debug "Launching marker build for $marker.codon\b";
		$job_id = launch_marker_build( marker => $marker, dna => 1, batch_script => $bs );
		$job_ids{$job_id} = 1 if defined $job_id;
	}
	my @taxonids = keys(%alltaxa);
	Phylosift::MarkerBuild::make_ncbi_subtree( out_file => "$marker_dir/ncbi_tree.updated.tre", taxon_ids => \@taxonids );

	wait_for_jobs( job_ids => \%job_ids );
}

sub clean_representatives {
	my %args    = @_;
	my $file    = $args{infile} || miss("infile");
	my $outfile = $args{outfile} || miss("outfile");
	my $INALN   = ps_open($file);
	my $OUTALN  = ps_open(">$outfile");
	while ( my $line = <$INALN> ) {
		unless ( $line =~ /^>/ ) {
			my $first_lower;
			my $last_upper;

			# strip all chars outside the first and last uppercase character -- these are the boundaries of the homologous region
			for ( my $i = 0; $i < 2; $i++ ) {
				if ( $line =~ /[A-Z]/ ) {
					$line = substr( $line, $-[0] );
				}
				$line = reverse($line);
			}
			$line =~ s/-//g;             # remove gap chars
			$line =~ s/\.//g;            # remove gap chars
			$line =~ tr/[a-z]/[A-Z]/;    # make all chars uppercase for lastal (it uses lowercase as soft-masking)
			$line .= "\n";               # trailing newline was removed above
		}
		print $OUTALN $line;
	}
}

sub read_gene_ids {
	my %args = @_;
	my $file = $args{file} || miss("file");
	my $IDS  = ps_open($file);
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
	`cd $repository/markers ; rm taxdmp.zip ; taxit new_database `;
}

sub make_marker_taxon_map {
	my %args      = @_;
	my $self      = $args{self} || miss("self");
	my $markerdir = $args{marker_dir} || miss("marker_dir");

	my $AAIDS = ps_open( "$markerdir/".get_gene_id_file( dna => 0 ) );
	my $MARKERTAXONMAP = ps_open(">$markerdir/marker_taxon_map.updated.txt");

	while ( my $line = <$AAIDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		print $MARKERTAXONMAP "$uniqueid\t$taxon\n";
	}
	close $MARKERTAXONMAP;
}

sub make_ncbi_tree_from_update {
	my %args      = @_;
	my $self      = $args{self} || miss("self");
	my $markerdir = $args{marker_dir} || miss("marker_dir");

	my $AAIDS = ps_open( "$markerdir/".get_gene_id_file( dna => 0 ) );
	my @taxonids;

	while ( my $line = <$AAIDS> ) {
		chomp $line;
		my ( $marker, $taxon, $uniqueid ) = split( /\t/, $line );
		push( @taxonids, $taxon ) if $taxon =~ /^\d+$/;
	}
	Phylosift::MarkerBuild::make_ncbi_subtree( out_file => "$markerdir/ncbi_tree.updated.tre", taxon_ids => \@taxonids );
}

sub get_marker_name_base {
	my %args       = @_;
	my $name       = $args{marker} || miss("marker");
	my $dna        = $args{dna} || 0;
	my $updated    = $args{updated} || 0;
	my $pruned     = $args{pruned} || 0;
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

sub create_taxon_id_table {
	my %args      = @_;
	my $alignment = $args{alignment} || miss("alignment");
	my $output    = $args{output} || miss("output");
	my $alltaxa   = $args{alltaxa} || miss("alltaxa");

	#	print "OUTPUT :$output\n";
	my $OUTPUT = ps_open(">$output");
	my $INALN  = ps_open($alignment);
	while ( my $line = <$INALN> ) {
		if ( $line =~ /^>(.+)/ ) {
			my $header = $1;
			my $ncbi   = 0;
			if ( $header =~ /\.(\d+?)\.fasta/ ) {

				#			if ( $header =~ /^(\d+?)\./ ) {
				$ncbi = 1;
				my $taxon_id = $1;
				$ncbi = 0 unless defined $taxon_id;
				$ncbi = 0 if $taxon_id < 100;
				$ncbi = 0 if is_unclassified( taxon_id => $taxon_id );
				if ($ncbi) {
					print $OUTPUT "$header\t$taxon_id\n";
					$alltaxa->{$1} = 1;
				}
			}
			if ( $ncbi == 0 ) {
				my $name = $1 if $header =~ /^(.+?)\.\d+.+/;
				print $OUTPUT "$header\t$name\n";
			}
		}
	}
	close $OUTPUT;
}

sub is_unclassified {
	my %args = @_;
	my $taxon_id = $args{taxon_id} || miss("taxon_id");
	return 0 unless looks_like_number($taxon_id);
	my $parent = Phylosift::Summarize::read_ncbi_taxonomy_structure();
	my $t      = $taxon_id;

	#	print $taxon_id."\n";
	# need to remove anything that NCBI considers "unclassified", currently taxon node ID 12908.
	for ( ; defined($t) && $t != 1; $t = $parent->{$t}->[0] ) {
		last if $t == 12908;
	}
	return defined($t) && $t == 12908;
}

sub filter_short_and_unclassified_seqs_from_fasta {
	my %args         = @_;
	my $input_fasta  = $args{input_fasta} || miss("input_fasta");
	my $output_fasta = $args{output_fasta} || miss("output_fasta");
	my $min_pct      = $args{min_pct} || miss("min_pct");

	# create a pruned fasta
	my $FASTA       = ps_open($input_fasta);
	my $PRUNEDFASTA = ps_open( ">".$output_fasta );

	#	print "OUTPUTFASTA $output_fasta\n";
	my $curhdr = "";
	my $curseq = "";
	my %known;
	while ( my $line = <$FASTA> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			if ( length($curseq) > 0 ) {
				my $glen = $curseq =~ tr/-/-/;
				my $known_pct = 100 * ( length($curseq) - $glen ) / length($curseq);

				#				print $curhdr."\n";
				my $taxon_id = $1 if $curhdr =~ /\.(\d+?)\.fasta/;
				$taxon_id = $1 if !defined($taxon_id) && $curhdr =~ /^>(\S+)\.fasta/;

				#				my $taxon_id = $1 if $curhdr =~ />(\d+?)\./;
				#				print "TAXON : $taxon_id\n";
				print $PRUNEDFASTA "$curhdr\n$curseq\n"
				  if $known_pct > $min_pct
					  && ( ( !is_unclassified( taxon_id => $taxon_id ) && !defined( $known{$curhdr} ) || $Phylosift::Settings::keep_paralogs ) );
				$known{$curhdr} = 1;
				$curseq = "";
			}
			$curhdr = $line;
		} else {
			$curseq .= $line;
		}
	}
	if ( length($curseq) > 0 ) {
		my $glen = $curseq =~ tr/-/-/;
		my $known_pct = 100 * ( length($curseq) - $glen ) / length($curseq);
		my $taxon_id = $1 if $curhdr =~ /\.(\d+?)\.fasta/;
		$taxon_id = $1 if !defined($taxon_id) && $curhdr =~ /^>(\S+)\.fasta/;

		#		my $taxon_id = $1 if $curhdr =~ />(\d+?)\./;
		print $PRUNEDFASTA "$curhdr\n$curseq\n"
		  if $known_pct > $min_pct && ( ( !is_unclassified( taxon_id => $taxon_id ) && !defined( $known{$curhdr} ) || $Phylosift::Settings::keep_paralogs ) );
	}
}

sub filter_fasta {
	my %args         = @_;
	my $input_fasta  = $args{input_fasta} || miss("input_fasta");
	my $output_fasta = $args{output_fasta} || miss("output_fasta");
	my $keep_taxa    = $args{keep_taxa} || miss("keep_taxa");

	# create a pruned fasta
	my $FASTA       = ps_open($input_fasta);
	my $PRUNEDFASTA = ps_open( ">".$output_fasta );
	my $printing    = 0;
	while ( my $line = <$FASTA> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			$printing = defined( $keep_taxa->{$1} ) ? 1 : 0;
		}
		print $PRUNEDFASTA "$line\n" if $printing;
	}
}

sub qsub_job {
	my %args      = @_;
	my $script    = $args{script} || miss("script");
	my $saref     = $args{script_args};
	my $qsub_args = $args{qsub_args} || "";
	my $jobsref   = $args{job_ids};

	#	my $qsub_cmd = "qsub -q all.q -q eisen.q $qsub_args $script ";
	my $qsub_cmd = "qsub $qsub_args $script ";
	$qsub_cmd .= join( " ", @$saref ) if defined($saref);
	my $job = `$qsub_cmd`;
	$job =~ /Your job (\d+) /;
	$jobsref->{$1} = 1;
}

sub write_marker_build_script {
	my %args        = @_;
	my $destination = $args{destination} || miss("destination");
	my $bs          = $args{batch_script};
	my $BUILDSCRIPT = ps_open(">$bs");
	my $ps          = $FindBin::Bin;
	print $BUILDSCRIPT <<EOF;
#!/bin/sh
#\$ -cwd
#\$ -V
#\$ -S /bin/bash

export PATH="\$PATH:$ps"
$ps/phylosift build_marker -f --debug --alignment=\$1 --update-only --taxonmap=\$2 --reps_pd=\$3 \$4 --destination=$destination

EOF

}

sub launch_marker_build {
	my %args   = @_;
	my $dna    = $args{dna};
	my $marker = $args{marker} || miss("marker");
	my $bs     = $args{batch_script} || miss("batch_script");

	my %jobids;

	my $marker_fasta = get_fasta_filename( marker => $marker, updated => 1, dna => $dna );
	print STDERR "Couldnt find $marker_fasta\n" unless -e $marker_fasta;
	return unless -e $marker_fasta;
	return if $marker =~ /^PMPROK/ || $marker =~ /^DNGNGWU/;    # skip these since they all go in the concat

	my $qsub_args = $marker eq "concat" ? " -l mem_free=20G -v OMP_NUM_THREADS=6 " : "";
	my @marray = ( $marker_fasta, $marker_fasta.".taxon_ids", "0.01" );
	my $clean_reps = get_reps_filename( marker => $marker, updated => 1, clean => 1 );
	push( @marray, "--unaligned=$clean_reps" ) if $dna == 0;
	qsub_job( script => $bs, qsub_args => $qsub_args, job_ids => \%jobids, script_args => \@marray );
	my @jabbar = keys(%jobids);
	return $jabbar[0];
}

=head2 update_rna

RNA markers contain sequences for many taxa that aren't available in genome databases.
We want the updated markers to retain these original sequences.

=cut

sub update_rna {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $marker = $args{marker} || miss("marker");

	#	my $id_file     = $args{id_file} || miss("id_file");
	my $marker_dir = $args{marker_dir} || miss("marker_dir");

	my $base_alignment = Phylosift::Utilities::get_marker_aln_file( self => $self, marker => $marker );
	my $updated_alignment = get_fasta_filename( marker => $marker, updated => 1 );

	#	$AAID = ps_open( ">>" . $id_file );
	my $INALN  = ps_open($base_alignment);
	my $OUTALN = ps_open( ">>$marker_dir/".$updated_alignment );
	my $max_id = 0;
	while ( my $line = <$INALN> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			my $countstring = sprintf( "%010u", $max_id );
			my $ggid = $1;
			$ggid =~ s/\/\d+-\d+//g;    # remove trailing bioperl rubbish

			#			print $AAID "gg_$ggid\t$countstring\n";
			$line = ">gg_$ggid";
			$max_id++;
		}
		print $OUTALN "$line\n";
	}

	#	close $AAID;

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
	my ( $idtref, $mtiref ) = read_gene_ids( file => "$marker_dir/".get_gene_id_file( dna => 0 ) );
	my %id_to_taxon        = %$idtref;
	my %marker_taxon_to_id = %$mtiref;

	croak("$constraint_tree does not exist") unless -e $constraint_tree;

	# first convert the constraint tree leaf names to match the target
	my $tree = Bio::Phylo::IO->parse(
									  '-file'   => $constraint_tree,
									  '-format' => 'newick',
	)->first;
	print STDERR "Tree has ".scalar( @{ $tree->get_entities } )." nodes\n";

	# collect all the sequence IDs present in the target alignment
	my %target_ids;
	my $target_alignment         = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned );
	my $target_alignment_amended = get_fasta_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ).".amended";
	my $column_count             = 0;
	my $TALN                     = ps_open($target_alignment);
	my $dat                      = "";
	while ( my $line = <$TALN> ) {
		chomp $line;
		if ( $line =~ /^>(.+)/ ) {
			$target_ids{$1} = 1;
			$column_count = length($dat) if ( $column_count == 0 );
		} elsif ( $column_count == 0 ) {
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
			push( @missing_taxa, $rename ) unless defined( $target_ids{$rename} );
		} else {
			push( @missing_taxa, $name );
		}
		$node->set_name($rename);
	}
	my $renamed_tree = "$constraint_tree.renamed";
	my $RENAMETREE   = ps_open(">$renamed_tree");
	print $RENAMETREE unparse( '-phylo' => $tree, '-format' => 'newick' )."\n";
	close $RENAMETREE;

	#
	# create an amended alignment with any taxa missing from the target data
	`cp "$target_alignment" "$target_alignment_amended"`;
	my $AMENDALN = ps_open(">>$target_alignment_amended");
	foreach my $missing (@missing_taxa) {
		print $AMENDALN ">$missing\n".( "-" x $column_count )."\n";
	}
	close $AMENDALN;

	unless ($copy_tree) {

		#
		# enumerate the splits in the constraint tree
		my $constraint_splits = "$renamed_tree.splits";
		my $constraint_cl     = "print_splits $renamed_tree $constraint_splits";
		system($constraint_cl);

		#
		# now infer a target tree that respects protein tree constraints
		my $target_constrained_tree = get_fasttree_tre_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ).".constrained";
		my $target_constrained_log = get_fasttree_log_filename( marker => $target_marker, dna => 0, updated => 1, pruned => $target_pruned ).".constrained";

		my $data_type_args = "";
		$data_type_args = "-nt -gtr" unless ( Phylosift::Utilities::is_protein_marker( marker => $target_marker ) );

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
		my %job_ids = ();
		qsub_job( script => "/tmp/constrained_tre.sh", job_ids => \%job_ids, qsub_args => "-l mem_free=12G" );
		wait_for_jobs( job_ids => \%job_ids );

		# finally, make a pplacer package with the constrained tree
		my $taxit_cl =
		  "taxit create -a \"Aaron Darling\" -d \"topology-constrained marker $target_marker\" -l $target_marker -f $target_alignment_amended -t $target_constrained_tree -s $target_constrained_log -P $target_marker.constrained";
		system($taxit_cl);
		unlink($constraint_splits);
	} else {

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
	make_constrained_tree(
						   constraint_marker => "16s_reps_bac",
						   constraint_tree   => $tctree,
						   target_marker     => "concat",
						   marker_dir        => $marker_dir,
						   constraint_pruned => 0,
						   target_pruned     => 1,
						   copy_tree         => 1
	);

	# with a little luck we can now??? do what?
}

sub get_gene_id_file {
	my %args = @_;
	my $dna = $args{dna} || 0;
	return "gene_ids.codon.txt" if $dna;
	return "gene_ids.aa.txt";
}

sub package_markers {
	my %args            = @_;
	my $marker_dir      = $args{marker_directory} || miss("marker_directory");
	my $base_marker_dir = $args{base_marker_directory} || miss("base_marker_directory");
	my @markerlist      = Phylosift::Utilities::gather_markers();
	unshift( @markerlist, "concat" );

	# first collect the base markers
	foreach my $marker (@markerlist) {

		# copy the base marker into our directory
		`cp -r $base_marker_dir/$marker/ $marker_dir`;
		`cp -r $base_marker_dir/$marker.short/ $marker_dir` if -d "$base_marker_dir/$marker.short/";
	}

	# clean up some stuff that really shouldn't be packaged
	chdir($marker_dir);
	system("rm -f *.fasttree.log");
	system("rm -f *.updated*fasta");
	system("rm -f *.reps.clean");
	system("rm -f *.updated.reps");
	system("rm -f core.*");
	system("rm -f ps_build_marker.*");
	system("rm -rf */temp_ref");

	chdir( $marker_dir."/../" );
	system("pwd");
	my @timerval = localtime();
	my $datestr  = Phylosift::Utilities::get_date_YYYYMMDD;
	system("tar czf markers_$datestr.tgz markers");
	`rm -f markers.tgz`;
	`ln -s markers_$datestr.tgz markers.tgz`;
}
1;
