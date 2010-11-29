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
    

    #need to convert the "seed" alignment to stockholm format due to the --mapali requirement
#    my $fastaAli = new Bio::AlignIO( -file   => "$workingDir/markers/$marker.ali",
#				     -format => 'fasta');


#    my $stockAli = new Bio::AlignIO( -file   => ">$workingDir/Amph_temp/alignments/$marker.seed.stock",
#				     -format => 'stockholm',
#				     -flush  => 0); # go as fast as we can!
    
#    while(my $aln = $fastaAli->next_aln() ) { $stockAli->write_aln($aln) };
 
    
    
    `perl $workingDir/fasta2stockholm.pl $workingDir/markers/$marker.trimfinal > $alignDir/$marker.seed.stock`;

#using hmmer3
    if(!-e "$alignDir/$marker.stock.hmm"){
	`hmmbuild $alignDir/$marker.stock.hmm $alignDir/$marker.seed.stock`;
    }

#using hmmer2
#    if(!-e "$workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm"){
#	`/usr/bin/hmmbuild $workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm $workingDir/Amph_temp/alignments/$marker.seed.stock`;
#    }


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
    
#using hmmer3
    `hmmalign --outformat afa -o $alignDir/$marker.aln_hmmer3.fasta --mapali $alignDir/$marker.seed.stock $alignDir/$marker.stock.hmm $blastDir/$marker.candidate`;

#using hmmer2
#    `/usr/bin/hmmalign -o $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock --mapali $workingDir/Amph_temp/alignments/$marker.seed.stock $workingDir/Amph_temp/alignments/$marker.stock_hmmer2.hmm $workingDir/Amph_temp/Blast_run/$marker.candidate`;
    
#    my $stockAli = new Bio::AlignIO( -file   => "$workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock",-format => 'stockholm');
#    my %lengths= ();


#    open(fastaOUT,">$workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta")or die "Couldn't open $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta for writing\n";
#    while(my $aln = $stockAli->next_aln()){
#	print "length : ".$aln->length."\n";
#	foreach my $seq ($aln->each_seq()){
#	    print fastaOUT ">".$seq->id."\n";
#	    print fastaOUT $seq->seq."\n";
#	    print ">".$seq->id."\n";
#	    print $seq->seq."\n";
#	}
#    }

    #does not work for what we want
    #`perl $workingDir/stockholm2fasta.pl $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.stock > $workingDir/Amph_temp/alignments/$marker.aln_hmmer2.fasta`;

    
    #trim the alignment
    
    #find out all the indexes that have a . in the reference sequences
    my $originAli = new Bio::AlignIO(-file=>"$workingDir/markers/$marker.trimfinal", -format=>'fasta');
    my %referenceSeqs = ();
    $refAliLength = 0;
    while(my $aln = $originAli->next_aln()){
	foreach my $seq($aln->each_seq()){
	    #print $seq->id."\n";
	    $referenceSeqs{$seq->id}=1;
	    $refAliLength = length($seq->seq);
#	    print $seq->id."\t".$refAliLength."\t";
	}
#	$refAliLength=$aln->length;
    }
    my $count = 0;
    my $numSeq=0;
    my %insertHash=();
    my $aliIN = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta", -format=>'fasta');
    while(my $aln = $aliIN->next_aln()){
	foreach my $seq ($aln->each_seq()){
	    if(exists $referenceSeqs{$seq->id}){
		#split a sequence into individual characters to find all the '.'
		my $i =0;
		my @characters = split(//,$seq->seq());
		$numSeq++;
		$count=0;
		foreach my $char (@characters){
		    #a hash value of 1 means a . was found at that index in all reference sequences
		    #a hash value of 0 means in at least 1 reference sequence a non . was found at that index
		    if ($char eq '.'){
#			print $char." ";
#			print "'".$char."'\t$i\t";
			$count++;
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
#		print $count."\n";
	    }#else do nothing
	}
    }
#    print "Found ".$count/$numSeq." .\n";
    my $hmmer3Ali = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta",-format=>'fasta');
    open(aliOUT,">$alignDir/$marker.aln_hmmer3.trim")or die "Couldn't open $alignDir/$marker.aln_hmmer3.trim for writting\n";
    while(my $aln = $hmmer3Ali->next_aln()){
        foreach my $seq ($aln->each_seq()){
            my @characters=();
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
	    
		my $newSeq = "@characters";
		$newSeq =~ s/\s//g;
		print "$newSeq\n";
		$newSeq =~ s/[\.\*]/\-/g;
		print "$newSeq\n";
#		print $seq->id."\tnewSeq = $newSeq\n";
#		print "$marker\tOriginal length: $refAliLength \tbefore TRIM : $strlength\t newLength : ".length($newSeq)."\n";
#		exit;

		my $newIDs = $seq->id;
		$newIDs =~ s/[^\w\d]/_/g;
		
		print "old\t".$seq->id."\t new\t$newIDs\n";
		
		print aliOUT ">".$newIDs."\n".$newSeq."\n";
		
	    }
        }
    }
    close(aliOUT);
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

#    $pm->finish;
}

#$pm->wait_all_children;


#rewrite the marker.list for the markers that have hits
open(markersOUT,">$fileDir/markers.list")or die "Couldn't open $fileDir/markers.list\n";
if(scalar(@markers) >0){
    foreach my $marker (@markers){
        if($marker){
            print markersOUT $marker."\n";
        }
    }
}else{
    
}
close(markersOUT);
