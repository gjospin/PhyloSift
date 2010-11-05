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

my $usage = qq~
Usage : $0 <options> <marker.list>

~;

my $threadNum=1;
GetOptions("threaded=i" => \$threadNum,
    ) || die $usage;


die $usage unless ($ARGV[0]);

my $workingDir = getcwd;
my $markersFile = $ARGV[0];

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
`mkdir $workingDir/Amph_temp` unless (-e "$workingDir/Amph_temp");

#check if the Blast_run directory exists, if it doesn't create it.
`mkdir $workingDir/Amph_temp/alignments` unless (-e "$workingDir/Amph_temp/alignments");

foreach my $marker (@markers){

    if(!-e "$workingDir/Amph_temp/Blast_run/$marker.candidate"){
	print STDERR "Couldn't find $workingDir/Amph_temp/Blast_run/$marker.candidate";
	next;
    }

    #need to convert the "seed" alignment to stockholm format due to the --mapali requirement
    my $fastaAli = new Bio::AlignIO( -file   => "$workingDir/markers/$marker.ali",
				    -format => "fasta");
    my $stockAli = new Bio::AlignIO( -file   => ">$workingDir/Amph_temp/alignments/$marker.seed.stock",
				     -format => "stockholm");
#				     -flush  => 0); # go as fast as we can!
				
    while($seq = $fastaAli->next_aln) { $stockAli->write_aln($seq) };



    `hmmalign --outformat afa -o $workingDir/Amph_temp/alignments/$marker.aln --mapali $workingDir/Amph_temp/alignments/$marker.seed.stock $workingDir/markers/$marker.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;

    #`hmmalign --outformat afa -o $workingDir/Amph_temp/alignments/$marker.aln --mapali $workingDir/markers/$marker.ali $workingDir/markers/$marker.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;
    exit;
}
