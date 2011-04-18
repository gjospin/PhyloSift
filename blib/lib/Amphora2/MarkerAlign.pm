package Amphora2::MarkerAlign;

use Cwd;
use warnings;
use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::SeqIO;
use List::Util qw(min);
#use Parallel::ForkManager;


=head1 NAME

Amphora2::MarkerAlign - Subroutines to align reads to marker HMMs

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Run HMMalign for a list of families from .candidate files located in $workingDir/Amph_temp/Blast_run/


input : Filename containing the marker list


output : An alignment file for each marker listed in the input file

Option : -threaded = #    Runs Hmmalign using multiple processors.


Perhaps a little code snippet.

    use Amphora2::Amphora2;

    my $foo = Amphora2::Amphora2->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 MarkerAlign

=cut

sub MarkerAlign {

	@ARGV = @_;
	my $usage = qq{
	Usage : $0 <options> <marker.list> <readsFile>

	};

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

	# markers with less than this fraction aligning to the model will be discarded
	my $alnLengthCutoff = 0.4;
	my $minAlignedResidues = 20;

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
		print STDERR "Couldn't find $blastDir/$marker.candidate\n";
		next;
	    }if(-s "$blastDir/$marker.candidate" == 0){
		print STDERR "WARNING : the candidate file for $marker is empty\n";
		delete $markers[$index];
		next;
	    }

	    #converting the marker's reference alignments from Fasta to Stockholm (required by Hmmer3)
		Amphora2::Utilities::fasta2stockholm( "$workingDir/markers/$marker.trimfinal", "$alignDir/$marker.seed.stock" );

	    #build the Hmm for the marker using Hmmer3
	    if(!-e "$alignDir/$marker.stock.hmm"){
		`$Amphora2::Utilities::hmmbuild $alignDir/$marker.stock.hmm $alignDir/$marker.seed.stock`;
	    }

	    #Running hmmsearch for the candidates against the HMM profile for the marker
	    `$Amphora2::Utilities::hmmsearch -E 0.1 --tblout $alignDir/$marker.hmmsearch.tblout $alignDir/$marker.stock.hmm $blastDir/$marker.candidate > $alignDir/$marker.hmmsearch.out`;
	    my %hmmHits=();
	    my %hmmScores=();
	    open(tbloutIN,"$alignDir/$marker.hmmsearch.tblout");
	    open(tblOUT,">>$alignDir/test.tblout");
	    my $countHits = 0;
	    while(<tbloutIN>){
		chomp($_);
		if($_ =~ m/^(\S+)\s+-\s+(\S+)\s+-\s+(\S+)\s+(\S+)/){
		    $countHits++;
		    print tblOUT "HIT : $1\tEVAL : $3\n";
		    my $hitname = $1;
		    my $basehitname = $1;
		    my $hitscore = $4;
		    # in case we're using 6-frame translation
		    $basehitname =~ s/_[fr][123]$//g;
		    if(!defined($hmmScores{$basehitname}) || $hmmScores{$basehitname} < $hitscore ){
			$hmmScores{$basehitname}=$hitscore;
			$hmmHits{$basehitname}=$hitname;
			print STDERR "Marker $marker : Adding $basehitname for $hitname with score $hitscore\n";
		    }else{
			print STDERR "Skipping $hitname because it's low scoring\n";
		    }
		}
	    }
	    close(tblOUT);
	    close(tbloutIN);
	    
	    # added a check if the hmmsearch found hits to prevent the masking and aligning from failing
	    if($countHits==0){
		print STDERR "WARNING : The hmmsearch for $marker found 0 hits, removing marker from the list to process\n";
		delete $markers[$index];
		next;
	    }

	    #reading the Hmmsearch output file
	    #my $alnIO = Bio::AlignIO->new(-format => "fasta", -file=>">$alignDir/$marker.newaln.faa");

	#    my $hmmsearchOUT = new Bio::SearchIO->new(-format => 'hmmer', -file => "$alignDir/$marker.hmmsearch.out");
	#    while(my $result = $hmmsearchOUT->next_hit){
	#	while( my $hit=$result->next_hit){
	#	    while(my $hsp = $hit->next_hsp){
	#		my $aln = $hsp->get_aln;
	#		my $alnIO = Bio::AlignIO->new(-format => "fasta", -file =>">$alignDir/$marker.newaln.faa");
	#		$alnIO->write_aln($aln);
	#	    }

	#	}

	#    }
	    open(newCandidate,">$alignDir/$marker.newCandidate");
	    my $seqin = new Bio::SeqIO('-file'=>"$blastDir/$marker.candidate");
	    while(my $sequence = $seqin->next_seq){
		my $baseid = $sequence->id;
		$baseid =~ s/_[fr][123]$//g;
		if(exists $hmmHits{$baseid} && $hmmHits{$baseid} eq $sequence->id){
		    print newCandidate ">".$sequence->id."\n".$sequence->seq."\n";
		}else{
		    print STDERR "skipping baseid $baseid seqid ".$sequence->id."\n";
		    print STDERR $hmmHits{$baseid}." was the top hit\n" if(exists $hmmHits{$baseid});
		}
	    }
	    close(newCandidate);

	    #Align the hits to the reference alignment using Hmmer3
	    `$Amphora2::Utilities::hmmalign --outformat afa -o $alignDir/$marker.aln_hmmer3.fasta --mapali $alignDir/$marker.seed.stock $alignDir/$marker.stock.hmm $alignDir/$marker.newCandidate`;

	    #trimming the alignment
	    
	    #find out all the indexes that have a . in the reference sequences
	    my $originAli = new Bio::AlignIO(-file=>"$workingDir/markers/$marker.trimfinal", -format=>'fasta');
	    my %referenceSeqs = ();
	   
	    while(my $aln = $originAli->next_aln()){
		foreach my $seq($aln->each_seq()){
		    $referenceSeqs{$seq->id}=1;
		}
	    }
	    my %insertHash=();
	    #reading the 
	    my $aliIN = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta", -format=>'fasta');
	    my $masqseq;
	    while(my $aln = $aliIN->next_aln()){
		# mask out the columns that didn't align to the marker, we'll want to remove
		# them below
		$masqseq = "\0" x $aln->length;
		foreach my $seq ($aln->each_seq()){
		    if(exists $referenceSeqs{$seq->id}){
			my $curseq = $seq->seq();
			my $ch1 = "\1";
			$curseq =~ s/\./$ch1/g;
			$curseq =~ s/\w/\0/g;
			$curseq =~ s/[\-\*]/\0/g;
			$masqseq |= $curseq;
		    }#else do nothing
		}
	    }
		# figure out which columns have the first and last marker data
		# nonmarker columns contain a . and were masked above
		my $firstCol = length($masqseq);
		my $ch1 = "\1";
		$masqseq =~ s/^$ch1+//g;
		$firstCol -= length($masqseq);
		$masqseq =~ s/$ch1//g;
		my $collen = length($masqseq);
	    # reading and trimming out non-marker alignment columns from Hmmalign output (Hmmer3)
	    my $hmmer3Ali = new Bio::AlignIO(-file =>"$alignDir/$marker.aln_hmmer3.fasta",-format=>'fasta');
	    open(aliOUT,">$alignDir/$marker.aln_hmmer3.trim")or die "Couldn't open $alignDir/$marker.aln_hmmer3.trim for writting\n";
	    while(my $aln = $hmmer3Ali->next_aln()){
		foreach my $seq ($aln->each_seq()){
		    #initialize empty array
		    my @characters=();
		    #if the sequence isn't a reference sequence remove all the characters at the indices that have a 1 value in the %insertHash
		    if(!exists $referenceSeqs{$seq->id}){
			# strip out all the columns that had a . (indicated by masqseq)
			my $newSeq = $seq->seq();
			my $seqLen= 0;
			#sequence length before masking
			$seqLen++ while $newSeq =~ m/[^-\.]/g;
			$newSeq = substr($newSeq, $firstCol, $collen);
			#change the remaining . into - for pplacer to not complain
			$newSeq =~ s/\./-/g;
			my $alignCount=0;
			$alignCount++ while $newSeq =~ m/[^-]/g;
			my $minRatio = min ($alignCount / $collen),($alignCount / $seqLen);
			print STDERR $seq->id."\t$minRatio\t$alignCount\n";
			next if ($minRatio < $alnLengthCutoff && $alignCount < $minAlignedResidues);
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

}


=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-amphora2-amphora2 at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Amphora2-Amphora2>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Amphora2::Amphora2


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Amphora2-Amphora2>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Amphora2-Amphora2>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Amphora2-Amphora2>

=item * Search CPAN

L<http://search.cpan.org/dist/Amphora2-Amphora2/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of Amphora2::MarkerAlign.pm
