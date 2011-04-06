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
print STDERR "checking\n";
#reading the list of markers
open(markersIN,"$markersFile") or die "Couldn't open the markers file\n";
while(<markersIN>){
    chomp($_);
    next if(-z  "$alignDir/$_.aln_hmmer3.trim");
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
    #adding a printed statement to check on the progress
    print STDERR "Running Placer on $marker ....\t";
    #running Pplacer
    if(!-e "$treeDir/$marker.aln_hmmer3.trim.place"){
	`pplacer -p -r $alignDir/$marker.trimfinal.fasta -t $workingDir/markers/$marker.final.tre -s $workingDir/markers/$marker.in_phyml_stats.txt $alignDir/$marker.aln_hmmer3.trim.fasta`;
    }
    #adding a printed statement to check on the progress (not really working if using parrallel jobs)
    print STDERR "Done !\n";
    #Pplacer write its output to the directory it is called from. Need to move the output to the trees directory
    if(-e "$workingDir/$marker.aln_hmmer3.trim.place"){
	`mv $workingDir/$marker.aln_hmmer3.trim.place $treeDir`;
    }

    #Transform the .place file into a tree file
    if(-e "$treeDir/$marker.aln_hmmer3.trim.place" && !-e "$treeDir/$marker.aln_hmmer3.trim.tog.tree"){
	`placeviz --loc -p $treeDir/$marker.aln_hmmer3.trim.place`
    }


    #placeviz writes its output to the directory it was called from, need to move the output to the trees directory
    my @placevizFiles = <$workingDir/$marker.aln_hmmer3.*>;
    if(scalar(@placevizFiles)>0){
	`mv $workingDir/$marker.* $treeDir`;
    }
#    if(-e "$workingDir/$marker.aln_hmmer3.*"){
#	`mv $workingDir/$marker.* $treeDir`;
#    }elsif(-e "$workingDir/$marker.aln_hmmer3.trim.PP.tog.tre"){
#	`mv $workingDir/$marker.* $treeDir`;
#    }


    #added the .PP. check to accomodate for what pplacer names its files (tax branch or master branch)
    # transform the taxon names in the tree file
    if(-e "$treeDir/$marker.aln_hmmer3.trim.tog.tre"){
	nameTaxa("$treeDir/$marker.aln_hmmer3.trim.tog.tre");
    }elsif(-e "$treeDir/$marker.aln_hmmer.trim.PP.tog.tre"){
	nameTaxa("$treeDir/$marker.aln_hmmer3.trim.PP.tog.tre");
    }

    $pm->finish
}
$pm->wait_all_children;



sub nameTaxa {
    my $filename = shift;

    # read in the taxon name map
    my %namemap;
    open(NAMETABLE, "$workingDir/markers/name.table")or die "Couldn't open $workingDir/markers/name.table\n";
    while( my $line = <NAMETABLE> ){
	chomp $line;
	my @pair = split(/\t/, $line);
	$namemap{$pair[0]} = $pair[1];
    }
    
    # parse the tree file to get leaf node names
    # replace leaf node names with taxon labels
    open( TREEFILE, $filename );
    my @treedata = <TREEFILE>;
    close TREEFILE;
    open( TREEFILE, ">$filename" );
    foreach my $tree(@treedata){
	my @taxanames = split( /[\(\)\,]/, $tree );
	foreach my $taxon( @taxanames ){
	    next if $taxon =~ /^\:/;        # internal node, no taxon label
	    my @taxondata = split( /\:/, $taxon );
	    next unless @taxondata > 0;
	    if(defined($namemap{$taxondata[0]})){
		my $commonName = $namemap{$taxondata[0]}."-".$taxondata[0];
		$tree =~ s/$taxondata[0]/$commonName/g;
	    }
	}
	print TREEFILE $tree;
    }
}

