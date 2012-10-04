#!/usr/bin/env perl
use strict;
use warnings;

# usage: prep_nextera.pl <input directory> <barcodes> <output directory>
#
# runs adapter trimming, screens for any vector/host contamination
# creates a pooled file in interleaved fasta for idba_ud, and also a set of interleaved fastq, one per sample
#

#
# TODO: this needs to get on our PhyloSift coding convention
#
my $r1_regex = "_R1_";
my $r2_regex = "_R3_";
my $i1_regex = "_R2_";
my $i2_regex = "_R4_";

my $input_dir = $ARGV[0];
my $barcode_file = $ARGV[1];
my $output_dir = $ARGV[2];

# the below need to be set appropriately
my $vector_db = "/home/koadman/data/test_next/vector";
my $adapter_db = "/home/koadman/data/test_next/adapter.fasta";

unless (-f $adapter_db){
	print STDERR "Error, adapter database $adapter_db not found. Please check the path and edit the script if needed\n";
	exit;
}

unless (-f "$vector_db.1.bt2"){
	print STDERR "Error, vector database $vector_db not found. Please check the path and edit the script if needed\n";
	exit;
}

my $min_reads = 1000;

my %r1files;
my %r2files;
my %bc_count;
read_barcodes(barcode_file => $barcode_file, r1files => \%r1files, r2files => \%r2files, output_dir => $output_dir);

my %r_files;	# stores input file names
get_input_files(input_dir => $input_dir, r_files=>\%r_files);
demultiplex(r1files=>\%r1files, r2files=>\%r2files, bc_count=>\%bc_count);
remove_unused_barcodes(r1files=>\%r1files, r2files=>\%r2files, bc_count=>\%bc_count);
vector_screen(r1files=>\%r1files, r2files=>\%r2files);
adapter_screen(r1files=>\%r1files, r2files=>\%r2files);

#
# all done!
#

sub read_fq_entry {
	my %args = @_;
	my $FILE = $args{FILE};
	return unless defined $FILE;
	return unless fileno $FILE;
	my @lines;
	for(my $i=0; $i<4; $i++){
		$lines[$i] = <$FILE>;
	}
	return @lines;
}

# organize the input file names in a datastructure
sub get_input_files {
	my %args = @_;
	my $input_dir = $args{input_dir};
	my $r_files = $args{r_files};
	open(LSSER, "ls -1 $input_dir |");
	while(my $line = <LSSER>){
		chomp $line;
		next unless $line =~ /fastq.gz/;
		if($line =~ /$r1_regex(\d+).fastq.gz/){
			$r_files->{$r1_regex}{$1} = $line;
		}elsif($line =~ /$r2_regex(\d+).fastq.gz/){
			$r_files->{$r2_regex}{$1} = $line;
		}elsif($line =~ /$i1_regex(\d+).fastq.gz/){
			$r_files->{$i1_regex}{$1} = $line;
		}elsif($line =~ /$i2_regex(\d+).fastq.gz/){
			$r_files->{$i2_regex}{$1} = $line;
		}
	}
}


sub read_barcodes {
	my %args = @_;
	my $barcode_file = $args{barcode_file};
	my $r1files = $args{r1files};
	my $r2files = $args{r2files};
	my $output_dir = $args{output_dir};

	open(ADAPTER, $barcode_file);
	my $bc_length = -1;
	while( my $line = <ADAPTER> ){
		chomp $line;
		my @adpt = split(/\t/, $line);
		my $cur_bc = $adpt[1];
		$bc_length = length($cur_bc) if $bc_length<0;
		my $cur_name = $adpt[0];
		unless( defined($r1files->{$cur_bc})){
			$r1files->{$cur_bc}{name} = "$output_dir/$cur_name.r1.fastq";
			$r2files->{$cur_bc}{name} = "$output_dir/$cur_name.r2.fastq";
			open($r1files->{$cur_bc}{FILE}, ">$r1files->{$cur_bc}{name}");
			open($r2files->{$cur_bc}{FILE}, ">$r2files->{$cur_bc}{name}");
		}
	}
}

sub demultiplex {
	my %args = @_;
#	my $r1files = $args{r1files};
#	my $r2files = $args{r2files};
#	my $bc_count = $args{bc_count};
	
	# read each input file and create de-multiplexed, interleaved fastq
	print STDERR "Demultiplexing barcodes\n";
	foreach my $chunk( keys(%{$r_files{$r1_regex}})){	
		open(my $R1, "zcat $input_dir/".$r_files{$r1_regex}{$chunk}." |") if defined($r_files{$r1_regex}{$chunk});
		open(my $R2, "zcat $input_dir/".$r_files{$r2_regex}{$chunk}." |") if defined($r_files{$r2_regex}{$chunk});
		open(my $I1, "zcat $input_dir/".$r_files{$i1_regex}{$chunk}." |") if defined($r_files{$i1_regex}{$chunk});
		open(my $I2, "zcat $input_dir/".$r_files{$i2_regex}{$chunk}." |") if defined($r_files{$i2_regex}{$chunk});
		my $counter = 0;
		while(1){
			my @i1_read = read_fq_entry(FILE=>$I1);
			my @i2_read = read_fq_entry(FILE=>$I2);
			my @r1_read = read_fq_entry(FILE=>$R1);
			my @r2_read = read_fq_entry(FILE=>$R2);
			
			last unless @i1_read && defined($i1_read[1]);
			my $barcode = $i1_read[1];
			chomp $barcode;
			next unless defined($r1files{$barcode});
	
			my $R1FILE = $r1files{$barcode}{FILE};
			my $R2FILE = $r2files{$barcode}{FILE};
			print $R1FILE @r1_read;
			print $R2FILE @r2_read;
			$bc_count{$barcode} = 0 unless defined $bc_count{$barcode};
			$bc_count{$barcode}++;
			$counter++;
		}
		print "Read $counter reads\n";
	}
	# close files, flush I/O buffers
	foreach my $barcode(keys(%r1files)){
			my $R1FILE = $r1files{$barcode}{FILE};
			my $R2FILE = $r2files{$barcode}{FILE};
			close $R1FILE;
			close $R2FILE;
	}
	
}

sub remove_unused_barcodes {
	my %args = @_;
	my $r1files = $args{r1files};
	my $r2files = $args{r2files};
	my $bc_count = $args{bc_count};
	# delete files with too few reads
	foreach my $barcode(keys(%$r1files)){
		next if defined($bc_count->{$barcode}) && $bc_count->{$barcode} > $min_reads;
		`rm -f $r1files->{$barcode}{name} $r2files->{$barcode}{name}`;
		delete $r1files->{$barcode};
		delete $r2files->{$barcode};
	}
}

sub vector_screen {
	my %args = @_;
	# screen each fastq for vector contamination
	print STDERR "Filtering out vector contamination\n";
	foreach my $barcode(keys(%r1files)){
		my $bwt_cmd = "bowtie2-align --local -p 8 -x $vector_db -1 ".$r1files{$barcode}{name}." -2 ".$r2files{$barcode}{name}." |";
		print STDERR "Running $bwt_cmd\n";
		open(my $BWT, $bwt_cmd);
		my $filt_name = $r1files{$barcode}{name};
		$filt_name =~ s/\.r1\.fastq/.filtered.fastq/;
		open(my $FILTFQ, ">$filt_name");
		my $pair_id=0;
		while(my $line = <$BWT>){
			next if $line =~ /^@/;
			chomp $line;
			my @vals = split(/\t/, $line);
			next unless $vals[2] eq "*";
			print $FILTFQ "\@$vals[0]/".($pair_id+1)."\n$vals[9]\n+\n$vals[10]\n";
			$pair_id++;
			$pair_id = $pair_id % 2;
		}
	}
}

##
## prep_nextera.pl /share/eisen-d5/illumina_data/pool_1a/slims.bioinformatics.ucdavis.edu/Data/5hy955z5db/Unaligned_3reads/Sample_Pool1_A/ /share/eisen-d5/illumina_data/plaque_miseq5/nextera_bcs.txt.csv /state/partition1/koadman/pool_1a/
## 
sub adapter_screen {
	# tagdust 'em?
	print STDERR "Filtering out adapter contamination\n";
	open(my $ALLCLEAN, ">$output_dir/all_pooled.fa");
	foreach my $barcode(keys(%r1files)){
		my $filt_name = $r1files{$barcode}{name};
		$filt_name =~ s/\.r1\.fastq/.filtered.fastq/;
		next unless -e $filt_name;
		my $td_name = $r1files{$barcode}{name};
		$td_name =~ s/\.r1\.fastq/.tagdust.fastq/;
		`mkfifo td.fifo`;
		`mkfifo sga.fifo`;
		my $tagdust_cmd = "tagdust -s -o td.fifo $adapter_db $filt_name & sga preprocess -q 15 -m 30 -o sga.fifo td.fifo &";
		system($tagdust_cmd);
		open(my $TD, "sga.fifo");
		open(my $TDOUT, ">$td_name");
		my @prev_read;
		while(1){
			my @read = read_fq_entry(FILE=>$TD);
			last unless @read && defined($read[1]);
			if(defined($prev_read[0])&&$prev_read[0] eq $read[0]){
				print $TDOUT @prev_read;
				print $TDOUT @read;
				my $fa_name = $prev_read[0];
				$fa_name =~ s/^@/>/g;
				chomp $fa_name;
				print $ALLCLEAN "$fa_name/1\n$prev_read[1]$fa_name/2\n$read[1]";
				@prev_read = ();
			}
			@prev_read = @read;
		}
		`rm -f td.fifo`;
		`rm -f sga.fifo`;
	}
}
	