#!/usr/bin/perl -w

#Run blast on a list of families 
#
# input : Filename with marker list
#         Filename for the reads file
#
#
# Output : For each marker, create a fasta file with all the reads and reference sequences.
#          Storing the files in a directory called Blast_run
#
# Option : -clean removes the temporary files created to run the blast
#          -threaded = #    Runs blast on multiple processors
#




use Getopt::Long;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;

my $usage = qq~
Usage: $0 <options> <marker.list> <reads_file>

~;

my $clean = 0;
my $threadNum = 1;

GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
    ) || die $usage;


die $usage unless ($ARGV[0] && $ARGV[1]);

my $workingDir = getcwd;
my $readsFile = $ARGV[1];
my $markersFile = $ARGV[0];

$readsFile =~ m/(\w+)\.(\w+)$/;
my $readsCore = $1;

my @markers = ();
#reading the list of markers
open(markersIN,"$markersFile") or die "Couldn't open the markers file\n";
while(<markersIN>){
    chomp($_);
    push(@markers, $_);
}
close(markersIN);

#check the input file exists
print "$readsFile was not found \n" unless (-e "$workingDir/$readsFile" || -e "$readsFile");

#check if the temporary directory exists, if it doesn't create it.
`mkdir $workingDir/Amph_temp` unless (-e "$workingDir/Amph_temp");

#check if the Blast_run directory exists, if it doesn't create it.
`mkdir $workingDir/Amph_temp/Blast_run` unless (-e "$workingDir/Amph_temp/Blast_run");

#remove rep.faa if it already exists
if(-e "$workingDir/Amph_temp/Blast_run/rep.faa"){ `rm $workingDir/Amph_temp/Blast_run/rep.faa`;}



#change the sequences to include the marker name as part of the ID
#also makes 1 large file with all the marker sequences
open(repOUT,">$workingDir/Amph_temp/Blast_run/rep.faa");
foreach my $marker (@markers){
    open(repIN,"$workingDir/markers/$marker.faa") or die "Couldn't open $marker.faa\n";
    while(<repIN>){
	chomp($_);
	if($_ =~ /^>(\S+)(\s.*)/){
	    print repOUT ">$1_$marker $2\n";
	    #	exit;
	}
	else{
	    print repOUT $_."\n";
	}
    }
    close(repIN);
}
close(repOUT);


get_blast_hits();
exit;

#make a blastable DB
`makeblastDB -in $workingDir/Amph_temp/Blast_run/rep.faa -dbtype prot -title RepDB`;


#print "Processing $coreReadsName\n";
#blast the reads to the DB
`blastp -query $readsFile -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $workingDir/Amph_temp/Blast_run/rep.faa -out $workingDir/amph_temp/Blast_run/$readsCore.blastp -outfmt 0 -num_threads $threadNum`;

my %markerHits = ();
get_blast_hits();

sub get_blast_hits{
    #parsing the blast file
    my %hits = ();
    my $in = new Bio::SearchIO('-format'=>'blast','-file' => "$workingDir/Amph_temp/Blast_run/$readsCore.blastp");
    while (my $result = $in->next_result) {
	#hit name is a markerName
	#print "Getting HIT \n";
	my $topFamily="";
	my $topScore=0;
	while( my $hit = $result->next_hit ) {
	    print $result->query_name."\t".$hit->name()."\t".$hit->raw_score."\t";
	    $hit->name() =~ m/([^_]+)$/;
	    my $markerHit = $1;
	    if($topFamily ne $markerHit){
		#compare the previous value
		#only keep the top hit
		if($topScore < $hit->raw_score){
		    $topFamily = $markerHit;
		    $topScore = $hit->raw_score;
		}#else do nothing
	    }#else do nothing
	} 
	$markerHits{$topFamily}{$result->query_name}=1;
	
    }
    unless (%markerHits) {
	system("rm $workingDir/Amph_temp/Blast_run/$readsCore.blastp");
	exit(1);
    }
    my $seqin = new Bio::SeqIO('-file'=>"$workingDir/$readsFile");
    my $seqout = new Bio::SeqIO('-file'=>">$workingDir/Amph_temp/Blast_run/$readsFile.candidate",'-format'=>'fasta');
   # while (my $seq = $seqin->next_seq) {
#		if (${$seq->id}) {
#			$seq{$seq->id} = $seq;
#			$seqout->write_seq($seq);
#		}
#	}

  
}
