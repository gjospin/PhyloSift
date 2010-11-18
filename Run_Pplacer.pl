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
Usage : $0 <options> <marker.list>

~;

#variable kept and set to 1 only until the multi thread is implemented.
my $threadNum = 1;

GetOptions("threaded=i" => \$threadNum)|| die $usage;


die $usage unless ($ARGV[0]);


my $workingDir = getcwd;
my $markersFile = $ARGV[0];

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
`mkdir $workingDir/Amph_temp` unless (-e "$workingDir/Amph_temp");

#check if the Blast_run directory exists, if it doesn't create it.
`mkdir $workingDir/Amph_temp/trees` unless (-e "$workingDir/Amph_temp/trees/");

my $pm = new Parallel::ForkManager($threadNum);

foreach my $marker(@markers){
    $pm->start and next;
    # Pplacer requires the alignment files to have a .fasta extension
    if(!-e "$workingDir/Amph_temp/alignments/$marker.trimfinal.fasta"){
	`cp $workingDir/markers/$marker.trimfinal $workingDir/Amph_temp/alignments/$marker.trimfinal.fasta`;
    }
    if(!-e "$workingDir/Amph_temp/alignments/$marker.aln_hmmer3.trim.fasta"){
        `cp $workingDir/Amph_temp/alignments/$marker.aln_hmmer3.trim $workingDir/Amph_temp/alignments/$marker.aln_hmmer3.trim.fasta`;
    }
    print STDERR "Running Placer on $marker ....\t";
    if(!-e "$workingDir/Amph_temp/trees/$marker.aln_hmmer3.trim.place"){
	`pplacer -p -r $workingDir/Amph_temp/alignments/$marker.trimfinal.fasta -t $workingDir/markers/$marker.final.tre -s $workingDir/markers/$marker.in_phyml_stats.txt $workingDir/Amph_temp/alignments/$marker.aln_hmmer3.trim.fasta`;
    }
    print STDERR "Done !\n";
    if(-e "$workingDir/$marker.aln_hmmer3.trim.place"){
	`mv $workingDir/$marker.aln_hmmer3.trim.place $workingDir/Amph_temp/trees`;
    }
    $pm->finish;
}
$pm->wait_all_children;
