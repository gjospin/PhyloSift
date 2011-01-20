#!/usr/bin/perl -w

#Run blast on a list of families for a set of Reads
#
# input : Filename with marker list
#         Filename for the reads file
#
#
# Output : For each marker, create a fasta file with all the reads and reference sequences.
#          Storing the files in a directory called Blast_run
#
# Option : -clean removes the temporary files created to run the blast
#          -threaded = #    Runs blast on multiple processors (I haven't see this use more than 1 processor even when specifying more)
#




use Getopt::Long;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;

my $usage = qq~
Usage: $0 <options> <marker.list> <reads_file>

~;

my $clean = 0; #option set up, but not used for later
my $threadNum = 1; #default value runs on 1 processor only.

GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
    ) || die $usage;


die $usage unless ($ARGV[0] && $ARGV[1]);

#set up filenames and directory variables
my $readsFile = $ARGV[1];
my $markersFile = $ARGV[0];
my %markerHits = ();

my $position = rindex($readsFile,"/");
my $fileName = substr($readsFile,$position+1,length($readsFile)-$position-1);

my $workingDir = getcwd;
my $tempDir = "$workingDir/Amph_temp";
my $fileDir = "$tempDir/$fileName";
my $blastDir = "$fileDir/Blast_run";

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
die "$readsFile was not found \n" unless (-e "$workingDir/$readsFile" || -e "$readsFile");

#check if the temporary directory exists, if it doesn't create it.
`mkdir $tempDir` unless (-e "$tempDir");

#create a directory for the Reads file being processed.
`mkdir $fileDir` unless (-e "$fileDir");

#check if the Blast_run directory exists, if it doesn't create it.
`mkdir $blastDir` unless (-e "$blastDir");

#remove rep.faa if it already exists (starts clean is a previous job was stopped or crashed or included different markers)
if(-e "$blastDir/rep.faa"){ `rm $blastDir/rep.faa`;}



#also makes 1 large file with all the marker sequences

foreach my $marker (@markers){
    #if a marker candidate file exists remove it, it is from a previous run and could not be related
    if(-e "$blastDir/$marker.candidate"){
	`rm $blastDir/$marker.candidate`;
    }

    #initiate the hash table for all markers incase 1 marker doesn't have a single hit, it'll still be in the results 
    #and will yield an empty candidate file
    $markerHits{$marker}="";
    
    #append the rep sequences for all the markers included in the study to the rep.faa file
    `cat $workingDir/markers/representatives/$marker.rep >> $blastDir/rep.faa`

}


#make a blastable DB
if(!-e "$blastDir/rep.faa.psq" ||  !-e "$blastDir/rep.faa.pin" || !-e "$blastDir/rep.faa.phr"){
    `makeblastdb -in $blastDir/rep.faa -dbtype prot -title RepDB`;
}


#check if a 6frame translation is needed
open(readCheck,$readsFile) or die "Couldn't open $readsFile\n";
my ($totalCount,$seqCount) =0;
while(<readCheck>){
    chomp($_);
    if($_=~m/^>/){
	next;
    }
    $seqCount++ while ($_ =~ /[atcgATCG]/g);
    $totalCount += length($_);
}
close(readCheck);
print STDERR "DNA % ".$seqCount/$totalCount."\n";
if($seqCount/$totalCount >0.8){
    print STDERR "Doing a six frame translation\n";
    #found DNA, translate in 6 frames
    `$workingDir/translateSixFrame $readsFile > $blastDir/$fileName-6frame`;
    if(!-e "$blastDir/$readsCore.blastp"){
	`blastp -query $blastDir/$fileName-6frame -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $blastDir/rep.faa -out $blastDir/$readsCore.blastp -outfmt 0 -num_threads $threadNum`;
    }
}else{
    #blast the reads to the DB
    if(!-e "$blastDir/$readsCore.blastp"){
	`blastp -query $readsFile -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $blastDir/rep.faa -out $blastDir/$readsCore.blastp -outfmt 0 -num_threads $threadNum`;
    }
}

get_blast_hits();

my %topscore = ();

sub get_blast_hits{
    #parsing the blast file
    my %hits = ();
    my $in = new Bio::SearchIO('-format'=>'blast','-file' => "$blastDir/$readsCore.blastp");
    while (my $result = $in->next_result) {
	#hit name is a markerName
	#print "Getting HIT \n";

	#keep track of the top family and score for each query
	my $topFamily="";
	my $topScore=0;
	while( my $hit = $result->next_hit ) {
	    #print $result->query_name."\t".$hit->name()."\t".$hit->raw_score."\t";
	    $hit->name() =~ m/([^_]+)$/;
	    my $markerHit = $1;
	    if($topFamily ne $markerHit){
		#compare the previous value
		#only keep the top hit
		if($topScore < $hit->raw_score){
		    $topFamily = $markerHit;
		    $topScore = $hit->raw_score;
		    $hits{$result->query_name}=$topFamily;
		}#else do nothing
	    }#else do nothing
	}
	$topscore{$result->query_name}=$topScore;
#	$markerHits{$topFamily}{$result->query_name}=1;

    }
    #if there are no hits, remove the blast file and exit (remnant from Martin's script)
    unless (%markerHits) {
	system("rm $blastDir/$readsCore.blastp");
	exit(1);
    }
    #read the readFile and grab the sequences for all the hit IDs and assign them to the right marker using a hash table
    my $seqin = new Bio::SeqIO('-file'=>"$readsFile");
    while (my $seq = $seqin->next_seq) {
	if(exists $hits{$seq->id}){
	    #print "\'".$hits{$seq->id}."\'\n";
	    #create a new string or append to an existing string for each marker
	    if(exists  $markerHits{$hits{$seq->id}}){
		$markerHits{$hits{$seq->id}} .= ">".$seq->id."_$topscore{$seq->id}\n".$seq->seq."\n";
	    }else{
		$markerHits{$hits{$seq->id}} = ">".$seq->id."_$topscore{$seq->id}\n".$seq->seq."\n";
	    }
	}
    }
    #write the read+ref_seqs for each markers in the list
    foreach my $marker (keys %markerHits){
	#writing the hits to the candidate file
	open(fileOUT,">$blastDir/$marker.candidate")or die "Couldn't open $blastDir/$marker.candidate for writing\n";
	print fileOUT $markerHits{$marker};
	close(fileOUT);
    }
}
