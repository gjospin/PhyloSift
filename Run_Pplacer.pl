#!/usr/bin/perl -w

#Runs Pplacer for all the markers in the list passed as input.
# the script will look for the files of interest in
# $workingDir/markers/
# AND
# $workingDir/Amph_temp/alignments/


use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Parallel::ForkManager;

my $usage = qq~
Usage : $0 <options> <marker.list> <ReadsFile>

~;

#variable kept and set to 1 only until the multi thread is implemented.
my $threadNum = 1;

GetOptions("threaded=i" => \$threadNum)|| die $usage;


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
my $treeDir = "$fileDir/trees";


my $markers=();

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
`mkdir $treeDir` unless (-e "$treeDir");

my $pm = new Parallel::ForkManager($threadNum);

foreach my $marker(@markers){
    $pm->start and next;
    # Pplacer requires the alignment files to have a .fasta extension
    if(!-e "$alignDir/$marker.trimfinal.fasta"){
	`cp $workingDir/markers/$marker.trimfinal $alignDir/$marker.trimfinal.fasta`;
    }
    if(!-e "$alignDir/$marker.aln_hmmer3.trim.fasta"){
        `cp $alignDir/$marker.aln_hmmer3.trim $alignDir/$marker.aln_hmmer3.trim.fasta`;
    }
    print STDERR "Running Placer on $marker ....\t";
    if(!-e "$treeDir/$marker.aln_hmmer3.trim.place"){
	`pplacer -p -r $alignDir/$marker.trimfinal.fasta -t $workingDir/markers/$marker.final.tre -s $workingDir/markers/$marker.in_phyml_stats.txt $alignDir/$marker.aln_hmmer3.trim.fasta`;
    }
    print STDERR "Done !\n";
    if(-e "$workingDir/$marker.aln_hmmer3.trim.place"){
	`mv $workingDir/$marker.aln_hmmer3.trim.place $workingDir/Amph_temp/$fileName/trees`;
    }
    $pm->finish
}
$pm->wait_all_children;
