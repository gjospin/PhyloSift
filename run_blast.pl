#!/usr/bin/perl

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


use warnings;
use strict;
use Getopt::Long;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;

my $usage = qq{
Usage: $0 <options> <marker.list> <reads_file>

};

my $usage2 = qq{
Usage : $0 <options> <marker.list> <reads_file_1> <reads_file_2>

};

my $clean = 0; #option set up, but not used for later
my $threadNum = 1; #default value runs on 1 processor only.
my $isolateMode=0; # set to 1 if running on an isolate assembly instead of raw reads
my $bestHitsBitScoreRange=30; # all hits with a bit score within this amount of the best will be used
my $pair=0; #used if using paired FastQ files
GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
	   "isolate" => \$isolateMode,
	   "bestHitBitScoreRange" => \$bestHitsBitScoreRange,
	   "paired" => \$pair,
    ) || die $usage;

my $readsFile_2 = "";
if($pair ==0){
    die $usage unless ($ARGV[0] && $ARGV[1]);
}elsif($pair != 0){
    die $usage2 unless ($ARGV[0] && $ARGV[1] && $ARGV[2]);
    $readsFile_2 = $ARGV[2];
}
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

$readsFile =~ m/(\w+)\.?(\w*)$/;
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
if($pair != 0){
    die "$readsFile_2 was not found \n" unless (-e "$workingDir/$readsFile_2" || -e "$readsFile_2");
}


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
#checking if paired fastQ files are used.
if($pair !=0){
    if(!-e "$blastDir/$fileName.fasta"){
	my %fastQ = ();
	my $curr_ID = "";
	my $skip = 0;
	print STDERR "Reading $readsFile\n";
	open(FASTQ_1, $readsFile)or die "Couldn't open $readsFile in run_blast.pl reading the FastQ file\n";
	while(<FASTQ_1>){
	    chomp($_);
	    if($_ =~ m/^@(\S+)\/(\d)/){
		$curr_ID =$1;
		$skip =0;
	    }elsif($_ =~ m/^\+$curr_ID\/(\d)/){
		$skip = 1;
	    }else{
		if($skip ==0){
		    $fastQ{$curr_ID}=$_;
		}else{
		    #do nothing
		}
	    }
	}
	close(FASTQ_1);
	print STDERR "Reading$readsFile_2 \n";
	open(FASTQ_2, $readsFile_2)or die "Couldn't open $readsFile_2 in run_blast.pl reading the FastQ file\n";
	while(<FASTQ_2>){
	    chomp($_);
	    if($_ =~ m/^@(\S+)\/(\d)/){
		$curr_ID =$1;
		$skip =0;
	    }elsif($_ =~ m/^\+$curr_ID\/(\d)/){
		$skip = 1;
	    }else{
		if($skip ==0){
		    my $reverse = reverse $_;
		    $fastQ{$curr_ID}.=$reverse;
		}else{
		    #do nothing
		}
	    }
	}
	close(FASTQ_2);
	print STDERR "Writing $fileName.fasta\n";
	open(FastA, ">$blastDir/$fileName.fasta")or die "Couldn't open $fileName.fasta for writing in run_blast.pl\n";
	foreach my $id (keys %fastQ){
	    print FastA ">".$id."\n".$fastQ{$id}."\n";
	}
	close(FastA);
    }
    #pointing $readsFile to the newly created fastA file
    $readsFile = "$blastDir/$fileName.fasta";
    if(!-e "$blastDir/$fileName-6frame"){
	`$workingDir/translateSixFrame $readsFile > $blastDir/$fileName-6frame`;
    }
    if(!-e "$blastDir/$readsCore.blastp"){
	`blastp -query $blastDir/$fileName-6frame -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $blastDir/rep.faa -out $blastDir/$readsCore.blastp -outfmt 6 -num_threads $threadNum`;
    }
    $readsFile = "$blastDir/$fileName-6frame";
}
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf STDERR "BEFORE 6Frame translation %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

if($pair ==0){
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
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR "Before 6Frame %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
	`$workingDir/translateSixFrame $readsFile > $blastDir/$fileName-6frame`;   
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	printf STDERR "After 6Frame %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
	if(!-e "$blastDir/$readsCore.blastp"){
	    `blastp -query $blastDir/$fileName-6frame -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $blastDir/rep.faa -out $blastDir/$readsCore.blastp -outfmt 6 -num_threads $threadNum`;
	}
	$readsFile="$blastDir/$fileName-6frame";
    }else{
    #blast the reads to the DB
	if(!-e "$blastDir/$readsCore.blastp"){
	    `blastp -query $readsFile -evalue 0.1 -num_descriptions 50000 -num_alignments 50000 -db $blastDir/rep.faa -out $blastDir/$readsCore.blastp -outfmt 6 -num_threads $threadNum`;
	}
    }
}
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf STDERR "After Blast %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
get_blast_hits();
my (%hitsStart,%hitsEnd, %topscore, %hits)=();
#my %topscore = ();
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf STDERR "After Blast Parse %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;


sub get_blast_hits{
    #parsing the blast file
    # parse once to get the top scores for each marker
    my %markerTopScores;
    my %topFamily=();
    my %topScore=();
    my %topStart=();
    my %topEnd=();
    open(blastIN,"$blastDir/$readsCore.blastp")or die "Couldn't open $blastDir/$readsCore.blastp\n";    
    while(<blastIN>){
	chomp($_);
	my @values = split(/\t/,$_);
	my $query = $values[0];
	my $subject = $values[1];
	my $query_start = $values[6];
	my $query_end = $values[7];
	my $bitScore = $values[11];
	my @marker = split(/\_/, $subject);
	my $markerName = $marker[$#marker];
	#parse once to get the top score for each marker (if isolate is ON, parse again to check the bitscore ranges)
	if($isolateMode==1){
	    # running on a genome assembly
	    # allow only 1 marker per sequence (TOP hit)
	    if( !defined($markerTopScores{$markerName}) || $markerTopScores{$markerName} < $bitScore ){
		$markerTopScores{$markerName} = $bitScore;
		$hitsStart{$query}{$markerName} = $query_start;
		$hitsEnd{$query}{$markerName}=$query_end;
	    }
#	}
#	if($isolateMode==1){
	    # running on a genome assembly
	    # allow more than one marker per sequence
	    # require all hits to the marker to have bit score within some range of the top hit
#	    if($markerTopScores{$markerName} < $hit->bits + $bestHitsBitScoreRange){
#		$hits{$hitName}{$markerHit}=1;
#		$hitsStart{$hitName}{$markerName} = $query_start;
#		$hitsEnd{$hitName}{$markerName} = $query_end;
#	    }
	}else{
	    # running on reads
	    # just do one marker per read
	    if(!exists $topFamily{$query}){
		$topFamily{$query}=$markerName;
		$topStart{$query}=$query_start;
		$topEnd{$query}=$query_end;
		$topScore{$query}=$bitScore;
	    }else{
		#only keep the top hit
		if($topScore{$query} <= $bitScore){
		    $topFamily{$query}= $markerName;
		    $topStart{$query}=$query_start;
		    $topEnd{$query}=$query_end;
		    $topScore{$query}=$bitScore;
		}#else do nothing
	    }#else do nothing
	}
    }
    close(blastIN);
    if($isolateMode ==1){
	# reading the output a second to check the bitscore ranges from the top score
	open(blastIN,"$blastDir/$readsCore.blastp")or die "Couldn't open $blastDir/$readsCore.blastp\n";
	# running on a genome assembly
	# allow more than one marker per sequence
	# require all hits to the marker to have bit score within some range of the top hit
	while(<blastIN>){
	    chomp($_);
	    my @values = split(/\t/,$_);
	    my $query = $values[0];
	    my $subject = $values[1];
	    my $query_start = $values[6];
	    my $query_end = $values[7];
	    my $bitScore = $values[11];
	    my @marker = split(/\_/, $subject);
	    my $markerName = $marker[$#marker];
	    if($markerTopScores{$markerName} < $bitScore + $bestHitsBitScoreRange){
		$hits{$query}{$markerName}=1;
		$hitsStart{$query}{$markerName} = $query_start;
		$hitsEnd{$query}{$markerName} = $query_end;
	    }
	}
	close(blastIN);
    }else{
	foreach my $queryID (keys %topScore){
	    $hits{$queryID}{$topScore{$queryID}}=1;
	    $hitsStart{$queryID}{$topScore{$queryID}}=$topStart{$queryID};
	    $hitsEnd{$queryID}{$topScore{$queryID}}=$topEnd{$queryID};
	}
    }

    my $seqin = new Bio::SeqIO('-file'=>"$readsFile");
    while (my $seq = $seqin->next_seq) {
	if(exists $hits{$seq->id}){
	    foreach my $markerHit(keys %{$hits{$seq->id}}){
		#print STDERR $seq->id."\t".$seq->description."\n";
		#checking if a 6frame translation was done and the suffix was appended to the description and not the sequence ID
		my $newID = $seq->id;
		if($seq->description =~ m/(_[fr][012])$/ && $seq->id !~m/(_[fr][012])$/){
		    $newID.=$1;
		}
		#create a new string or append to an existing string for each marker

		#pre-trimming for the query + 150 residues before and after (for very long queries)
		my $start = $hitsStart{$seq->id}{$markerHit}-150;
		if($start < 0){
		    $start=0;
		}
		my $end = $hitsEnd{$seq->id}{$markerHit}+150;
		my $seqLength = length($seq->seq);
		if($end >= $seqLength){
		    $end=$seqLength;
		}
		my $newSeq = substr($seq->seq,$start,$end-$start);
		if(exists  $markerHits{$markerHit}){
		    $markerHits{$markerHit} .= ">".$newID."\n".$newSeq."\n";
		}else{
		    $markerHits{$markerHit} = ">".$newID."\n".$newSeq."\n";
		}
	    }
	}
    }
    #write the read+ref_seqs for each markers in the list
    foreach my $marker (keys %markerHits){
	#writing the hits to the candidate file
	open(fileOUT,">$blastDir/$marker.candidate")or die " Couldn't open $blastDir/$marker.candidate for writing\n";
	print fileOUT $markerHits{$marker};
	close(fileOUT);
    }
}
