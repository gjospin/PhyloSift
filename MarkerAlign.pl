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

    if(!-e "$workingDir/Amph_temp/Blast_run/$marker.candidate" ){
	print STDERR "Couldn't find $workingDir/Amph_temp/Blast_run/$marker.candidate";
	next;
    }if(-s "$workingDir/Amph_temp/Blast_run/$marker.candidate" == 0){
	print STDERR "WARNING : the candidate file for $marker is empty\n";
	next;
    }
    

    #need to convert the "seed" alignment to stockholm format due to the --mapali requirement
#    my $fastaAli = new Bio::AlignIO( -file   => "$workingDir/markers/$marker.ali",
#				     -format => 'fasta');


#    my $stockAli = new Bio::AlignIO( -file   => ">$workingDir/Amph_temp/alignments/$marker.seed.stock",
#				     -format => 'stockholm',
#				     -flush  => 0); # go as fast as we can!
    
#    while(my $aln = $fastaAli->next_aln() ) { $stockAli->write_aln($aln) };
 
    
    
    `perl $workingDir/fasta2stockholm.pl $workingDir/markers/$marker.ali > $workingDir/Amph_temp/alignments/$marker.seed.stock`;

    `hmmbuild $workingDir/Amph_temp/alignments/$marker.stock.hmm $workingDir/Amph_temp/alignments/$marker.seed.stock`;

    if(!-e "$workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm"){
	`/usr/bin/hmmbuild $workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm $workingDir/Amph_temp/alignments/$marker.seed.stock`;
    }
#    open(stockTest, "$workingDir/Amph_temp/alignments/$marker.seed.stock")or die "Couldn't open the stockholm file for marker $marker";
#    open(stockOUT,">$workingDir/Amph_temp/alignments/$marker.seed.stock_2")or die "Couldn't open the stockholm output file";
#    while(<stockTest>){
#	chomp($_);
#	if($_ =~ m/^\#=GS/){
#	    $_ =~ s/\s{2,}//;
#	    $_ =~ s/unknown//;
#	    $_ =~ s/AC\s+/ AC /;
#	    print stockOUT $_."\n";
#	}else{
#	    print stockOUT $_."\n";
#	}
#    }
#    print stockOUT "//\n";
#    close(stockTest);
#    close(stockOUT);
    

    `hmmalign --outformat afa -o $workingDir/Amph_temp/alignments/$marker.aln_hmmer3.fasta --mapali $workingDir/Amph_temp/alignments/$marker.seed.stock $workingDir/Amph_temp/alignments/$marker.stock.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;

    `/usr/bin/hmmalign -o $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock --mapali $workingDir/Amph_temp/alignments/$marker.seed.stock $workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;
    
    my $stockAli = new Bio::AlignIO( -file   => "$workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock",-format => 'stockholm');
    my %lengths= ();


    open(fastaOUT,">$workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta")or die "Couldn't open $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta for writing\n";
    while(my $aln = $stockAli->next_aln()){
	print "length : ".$aln->length."\n";
	foreach my $seq ($aln->each_seq()){
	    print fastaOUT ">".$seq->id."\n";
	    print fastaOUT $seq->seq."\n";
#	    print ">".$seq->id."\n";
#	    print $seq->seq."\n";
	}
    }

    #does not work for what we want
    #`perl $workingDir/stockholm2fasta.pl $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock > $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta`;

    
    


    #check the file size, if no hits for the marker, the file size will be 0.
#    my $candidateSize = -s "$workingDir/Amph_temp/Blast_run/$marker.candidate";
#    print $marker."\t".$candidateSize."\n";
#    if($candidateSize != 0){
#	`hmmalign --outformat afa -o $workingDir/Amph_temp/alignments/$marker.hits_ali $workingDir/markers/$marker.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;

#	my $refAli = new Bio::AlignIO(-file => "$workingDir/markers/$marker.ali", -format => 'fasta');
#	my $hitAli = new Bio::AlignIO(-file => "$workingDir/Amph_temp/alignments/$marker.hits_ali", -format => 'fasta');
#	my $totalAli = new Bio::AlignIO(-file => ">$workingDir/Amph_temp/alignments/$marker.all_ali", -format =>'fasta');
#	while(my $aln = $refAli->next_aln() ) { $totalAli->write_aln($aln) };
#	while(my $aln = $hitAli->next_aln() ) { $totalAli->write_aln($aln) };
#    }
}


