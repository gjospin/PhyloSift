#!/usr/bin/perl -w

#Run HMMalign for a list of families from .candidate files located in $workingDir/Amph_temp/Blast_run/
#
#
# input : Filename containing the marker list
#
#
# output : An alignment file for each marker listed in the input file
#
# Option : -threaeded = #    Runs Hmmalign using multiple processors.
#
#
#
#
#



use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Parallel::ForkManager;


my $usage = qq~
Usage : $0 <options> <marker.list> <readsFile>

~;

my $threadNum=1;
GetOptions("threaded=i" => \$threadNum,
    ) || die $usage;


die $usage unless ($ARGV[0] && $ARGV[1]);

my $workingDir = getcwd;
my $markersFile = $ARGV[0];

my $readsFile = $ARGV[1];
my $position = rindex($readsFile,"/");
my $fileName = substr($readsFile,$position+1,length($readsFile)-$position-1);

my $tempDir = "$workingDir/Amph_temp";
my $fileDir = "$tempDir/$fileName";
my $blastDir = "$fileDir/Blast_run";
my $alignDir = "$fileDir/alignments";


my @markers = ();
#reading the list of markers
open(markersIN,"$markersFile") or die "Couldn't open the markers file\n";
while(<markersIN>){
    chomp($_);
    #print "\'$_\'\n";
    push(@markers, $_);
}
close(markersIN);

#check if the temporary directory exists, if it doesn't create it.
`mkdir $tempDir` unless (-e "$tempDir");

#create a directory for the Reads file being processed.
`mkdir $fileDir` unless (-e "$fileDir");

#check if the Blast_run directory exists, if it doesn't create it.
`mkdir $alignDir` unless (-e "$alignDir");
my $index =-1;

#my $pm = new Parallel::ForkManager($threadNum);
foreach my $marker (@markers){
#    $pm->start and next;
    $index++;
    if(!-e "$blastDir/$marker.candidate" ){
	print STDERR "Couldn't find $blastDir/$marker.candidate";
	next;
    }if(-s "$blastDir/$marker.candidate" == 0){
	print STDERR "WARNING : the candidate file for $marker is empty\n";
	delete $markers[$index];
	next;
    }
    
    #converting the marker's reference alignments from Fasta to Stockholm (required by Hmmer3)
    `perl $workingDir/fasta2stockholm.pl $workingDir/markers/$marker.trimfinal > $alignDir/$marker.seed.stock`;

    #build the Hmm for the marker using Hmmer3
    if(!-e "$alignDir/$marker.stock.hmm"){
	`hmmbuild $alignDir/$marker.stock.hmm $alignDir/$marker.seed.stock`;
    }

    #Align the hits to the reference alignment using Hmmer3
    `hmmalign --outformat afa -o $alignDir/$marker.aln_hmmer3.fasta --mapali $alignDir/$marker.seed.stock $alignDir/$marker.stock.hmm $blastDir/$marker.candidate`;

    #trimming the alignment
    
    #find out all the indexes that have a . in the reference sequences
    my $originAli = new Bio::AlignIO(-file=>"$workingDir/markers/$marker.trimfinal", -format=>'fasta');
    my %referenceSeqs = ();
    $refAliLength = 0;
    while(my $aln = $originAli->next_aln()){
	foreach my $seq($aln->each_seq()){
	    $referenceSeqs{$seq->id}=1;
	}
    }
    my %insertHash=();
    #reading the 
    my $aliIN = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta", -format=>'fasta');
    while(my $aln = $aliIN->next_aln()){
	foreach my $seq ($aln->each_seq()){
	    if(exists $referenceSeqs{$seq->id}){
		#split a sequence into individual characters to find all the '.'
		my $i =0;
		my @characters = split(//,$seq->seq());
		foreach my $char (@characters){
		    #a hash value of 1 means a . was found at that index in all reference sequences
		    #a hash value of 0 means in at least 1 reference sequence a non . was found at that index
		    if ($char eq '.'){
			if(!exists $insertHash{$i}){
			    $insertHash{$i}=1;
			}else{
			    $insertHash{$i} = $insertHash{$i}*1; #if all . have been found at that index, then the end will be 1
			}
		    }else{
			if(!exists $insertHash{$i}){
			    $insertHash{$i}=0;
			}else{
			    $insertHash{$i}= $insertHash{$i}*0;#if 1 non . character has been found at this index then it'll turn to 0
			}
		    }#else do nothing
		    $i++;
		}
	    }#else do nothing
	}
    }
    #reading and trimming the alignment from Hmmalign (Hmmer3)
    my $hmmer3Ali = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta",-format=>'fasta');
    open(aliOUT,">$alignDir/$marker.aln_hmmer3.trim")or die "Couldn't open $alignDir/$marker.aln_hmmer3.trim for writting\n";
    while(my $aln = $hmmer3Ali->next_aln()){
        foreach my $seq ($aln->each_seq()){
	    #initialize empty array
            my @characters=();
	    #if the sequence isn't a reference sequence remove all the characters at the indices that have a 1 value in the %insertHash
	    if(!exists $referenceSeqs{$seq->id}){
		my $i =0;
                @characters = split(//,$seq->seq());
		my $strlength = scalar(@characters);
#		print "STRlength : $strlength\t";
		while($i < $strlength){
		    if($insertHash{$i}==1){
			$characters[$i]='';
		    }
		    $i++;
                }
		#make a new sequence
		my $newSeq = "@characters";
		#remove the spaces introduced from the previous line
		$newSeq =~ s/\s//g;
		#transform . and * into -
		$newSeq =~ s/[\.\*]/\-/g;
		
		my $newIDs = $seq->id;
		#subsitute all the non letter or number characters into _ in the IDs to avoid parsing issues in tree viewing programs or others
		$newIDs =~ s/[^\w\d]/_/g;
		#print the new trimmed alignment
		print aliOUT ">".$newIDs."\n".$newSeq."\n";
		
	    }
        }
    }
    close(aliOUT);
#    $pm->finish;
}

#$pm->wait_all_children;


#rewrite the marker.list for the markers that have hits (the markers with no hits are removed from the list at this point)
open(markersOUT,">$fileDir/markers.list")or die "Couldn't open $fileDir/markers.list\n";
if(scalar(@markers) >0){
    foreach my $marker (@markers){
        if($marker){
            print markersOUT $marker."\n";
        }
    }
}else{
    #need to add a check with a warning if there are NO hits at all (the check in run_blast.pl shouldn't allow the script to get here anyways)
}
close(markersOUT);
