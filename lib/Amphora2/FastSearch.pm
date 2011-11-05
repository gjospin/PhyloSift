package Amphora2::FastSearch;

use warnings;
use strict;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::SeqUtils;
use Carp;
use Amphora2::Amphora2;
use Amphora2::Utilities qw(:all);
use File::Basename;
use POSIX qw(ceil floor);

use constant FLANKING_LENGTH => 150;

=head1 NAME

Amphora2::FastSearch - Subroutines to perform fast sequence identity searches between reads and marker genes.
Currently uses either BLAST or RAPsearch.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Run blast on a list of families for a set of Reads
 
 input : Filename with marker list
         Filename for the reads file


 Output : For each marker, create a fasta file with all the reads and reference sequences.
          Storing the files in a directory called Blast_run

 Option : -clean removes the temporary files created to run the blast
          -threaded = #    Runs blast on multiple processors (I haven't see this use more than 1 processor even when specifying more)

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 RunBlast

=cut

my $clean = 0; #option set up, but not used for later
my $threadNum = 4; #default value runs on 1 processor only.
my $isolateMode=0; # set to 1 if running on an isolate assembly instead of raw reads
my $bestHitsBitScoreRange=30; # all hits with a bit score within this amount of the best will be used
my $pair=0; #used if using paired FastQ files
my @markers;
my (%hitsStart,%hitsEnd, %topscore, %hits, %markerHits,%markerNuc)=();
my $readsCore;
my $custom="";
my %marker_lookup=();
my %frames=();
my $reverseTranslate=0;

my $blastdb_name = "blastrep.faa";
my $blastp_params = "-p blastp -e 0.1 -b 50000 -v 50000 -a $threadNum -m 8";

sub RunSearch {
    my $self = shift;
    my $custom = shift;
    my $searchtype = shift;
    my $markersRef = shift;
    @markers = @{$markersRef};
    %markerHits = ();
    my $position = rindex($self->{"readsFile"},"/");
    $self->{"readsFile"} =~ m/(\w+)\.?(\w*)$/;
    $readsCore = $1;
    $isolateMode = $self->{"isolate"};
    $reverseTranslate = $self->{"reverseTranslate"};
    debug "before rapPrepandclean\n";
    prepAndClean($self);
    if($self->{"readsFile_2"} ne ""){
	debug "before fastqtoFASTA\n";
	fastqToFasta($self);
    }

    my $resultsfile;
    if($searchtype eq "blast"){
	$resultsfile = blastXoof_table($self, $self->{"readsFile"});
	$reverseTranslate=1;
    }else{
        $resultsfile = executeRap($self);
	build_lookup_table($self);
    }
    get_hits($self,$resultsfile, $searchtype);
    return $self;
}

=head2 build_lookup_table

=cut

sub build_lookup_table{
    my $self=shift;
    debug "Building the lookup table for all markers";
    foreach my $markName (@markers){
	open(markIN,$Amphora2::Utilities::marker_dir."/".$markName.".faa");
	while(<markIN>){
	    chomp($_);
	    if ($_ =~ m/^>(\S+)/){
		$marker_lookup{$1}=$markName;
	    }
	}
	close(markIN);
	debug ".";
    }
    debug "\n";
    return $self;
}

=head2 blastXoof_table

=cut

sub blastXoof_table{
    my $self = shift;
    my $query_file = shift;
       debug "INSIDE tabular OOF blastx\n";
       `$Amphora2::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50000 -v 50000 -d $self->{"blastDir"}/$blastdb_name -o $self->{"blastDir"}/$readsCore.tabblastx -m 8 -a $threadNum`;
     return $self->{"blastDir"}."/$readsCore.tabblastx";
 }

=head2 blastXoof_full

=cut

sub blastXoof_full{
    my $query_file = shift;
    my $self = shift;
       debug "INSIDE full OOF blastx\n";
       `$Amphora2::Utilities::blastall -p blastx -i $query_file -e 0.1 -w 20 -b 50000 -v 50000 -d $self->{"blastDir"}/$blastdb_name -o $self->{"blastDir"}/$readsCore.blastx -a $threadNum`;
    return $self;
}

=head2 translateFrame

=cut

sub translateFrame{
    my $id = shift;
    my $seq = shift;
    my $start = shift;
    my $end = shift;
    my $frame = shift;
    my $marker = shift;
    my $reverseTranslate = shift;
    my $returnSeq = "";
    my $localseq = substr($seq, $start-1, $end-$start+1);
    
    my $newseq = Bio::LocatableSeq->new( -seq => $localseq, -id => 'temp');
    $newseq = $newseq->revcom() if($frame<0);
    if($reverseTranslate){
	$id =~ s/[\.\:\/\-]/_/g;
	if(exists  $markerNuc{$marker}){
	    $markerNuc{$marker} .= ">".$id."\n".$newseq->seq."\n";
	}else{
	    $markerNuc{$marker}= ">".$id."\n".$newseq->seq."\n";
	}
    }

    $returnSeq = $newseq->translate();

    return $returnSeq->seq();
}

=head2 executeRap

=cut

sub executeRap{
    my $self = shift;
    if($self->{"readsFile"} !~ m/^\//){
	debug "Making sure rapsearch can find the readsfile\n";
	$self->{"readsFile"}=getcwd()."/".$self->{"readsFile"};
	debug "New readsFile ".$self->{"readsFile"}."\n";
    }
    my $dbDir = "$Amphora2::Utilities::marker_dir/representatives";
    $dbDir = $self->{"blastDir"} if($custom ne "");
	if(!-e $self->{"blastDir"}."/$readsCore.rapSearch.m8"){
	    debug "INSIDE custom markers RAPSearch\n";
	    `cd $self->{"blastDir"} ; $Amphora2::Utilities::rapSearch -q $self->{"readsFile"} -d $dbDir/rep -o $readsCore.rapSearch -e -1`;
	}
    return $self->{"blastDir"}."/".$readsCore.".rapSearch.m8";
}

=head2 executeBlast

=cut

sub executeBlast{
    my $self = shift;
    my $query_file = shift;
    my $dbDir = "$Amphora2::Utilities::marker_dir/representatives";
    $dbDir = $self->{"blastDir"} if $custom ne "";
    debug "INSIDE BLAST\n";
    if(!-e $self->{"blastDir"}."/$readsCore.blastp"){
	`$Amphora2::Utilities::blastall $blastp_params -i $query_file -d $dbDir/$blastdb_name -o $self->{"blastDir"}/$readsCore.blastp`;
    }
    return $self->{"blastDir"}."/$readsCore.blastp";
}


=head2 fastqToFasta

    Writes a fastA file from 2 fastQ files from the Amphora2 object

=cut

sub fastqToFasta{
    my $self = shift;
    if($self->{"readsFile_2"} ne ""){
	debug "FILENAME ".$self->{"fileName"}."\n";

	return $self if(-e $self->{"blastDir"}."/$readsCore.fasta");
	
	my %fastQ = ();
	my $curr_ID = "";
	my $skip = 0;
	debug "Reading ".$self->{"readsFile"}."\n";
	open(FASTQ_1, $self->{"readsFile"})or die "Couldn't open ".$self->{"readsFile"}." in run_blast.pl reading the FastQ file\n";
	open(FASTQ_2, $self->{"readsFile_2"})or die "Couldn't open ".$self->{"readsFile_2"}." in run_blast.pl reading the FastQ file\n";            
	debug "Writing ".$readsCore.".fasta\n";
	open(FASTA, ">".$self->{"blastDir"}."/$readsCore.fasta")or die "Couldn't open ".$self->{"blastDir"}."/$readsCore.fasta for writing in run_blast.pl\n";
	while(my $head1 = <FASTQ_1>){
		my $read1 = <FASTQ_1>;
		my $qhead1 = <FASTQ_1>;
		my $qval1 = <FASTQ_1>;
		my $head2 = <FASTQ_2>;
		my $read2 = <FASTQ_2>;
		my $qhead2 = <FASTQ_2>;
		my $qval2 = <FASTQ_2>;

		$head1 =~ s/^\@/\>/g;
		chomp($read1);
		chomp($read2);
		$read2 =~ tr/ACGTacgt/TGCAtgca/;
		$read2 = reverse($read2);
		print FASTA "$head1$read1$read2\n";
	}

	#pointing $readsFile to the newly created fastA file
	$self->{"readsFile"} = $self->{"blastDir"}."/$readsCore.fasta";
    }
    return $self;
}

=head2 get_hits

parse the blast file

=cut



sub get_hits{
    my %duplicates = ();
    my $self = shift;
    my $hitfilename=shift;
    my $searchtype = shift;
    #parsing the blast file
    # parse once to get the top scores for each marker
    my %markerTopScores;
    my %topFamily=();
    my %topScore=();
    my %topStart=();
    my %topEnd=();
    open(blastIN,$hitfilename)or carp("Couldn't open ".$hitfilename."\n");
    while(<blastIN>){
	chomp($_);
	next if($_ =~ /^#/);
	my @values = split(/\t/,$_);
	my $query = $values[0];
	my $subject = $values[1];
	my $query_start = $values[6];
	my $query_end = $values[7];
	my $bitScore = $values[11];
	my $markerName = "";
	if($searchtype eq "blast"){
		my @marker=split(/\_/,$subject);
		$markerName = $marker[$#marker];
	}else{
		my @marker=split(/\_\_/,$subject);
		$markerName = $marker[0];
		
#		$markerName = $marker_lookup{$subject};
	}

	#parse once to get the top score for each marker (if isolate is ON, parse again to check the bitscore ranges)
	if($isolateMode==1){
	    # running on a genome assembly
	    # allow only 1 marker per sequence (TOP hit)
	    if( !defined($markerTopScores{$markerName}) || $markerTopScores{$markerName} < $bitScore ){
		$markerTopScores{$markerName} = $bitScore;
		$hitsStart{$query}{$markerName} = $query_start;
		$hitsEnd{$query}{$markerName}=$query_end;
	    }
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
	open(blastIN,$hitfilename)or die "Couldn't open $hitfilename\n";
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
		my $markerName = "";
		if($searchtype eq "blast"){
			my @marker=split(/\_/,$subject);
			$markerName = $marker[$#marker];
		}else{
			$markerName = $marker_lookup{$subject};
		}
	    if($markerTopScores{$markerName} < $bitScore + $bestHitsBitScoreRange){
		$hits{$query}{$markerName}=1;
		$hitsStart{$query}{$markerName} = $query_start;
		$hitsEnd{$query}{$markerName} = $query_end;
	    }
	}
	close(blastIN);
    }else{
	foreach my $queryID (keys %topFamily){
	    $hits{$queryID}{$topFamily{$queryID}}=1;
	    $hitsStart{$queryID}{$topFamily{$queryID}}=$topStart{$queryID};
	    $hitsEnd{$queryID}{$topFamily{$queryID}}=$topEnd{$queryID};
	}
    }
    my $seqin;
    if(-e $self->{"blastDir"}."/$readsCore-6frame"){
	$seqin = new Bio::SeqIO('-file'=>$self->{"blastDir"}."/$readsCore-6frame");
    }else{
	debug "ReadsFile:  $self->{\"readsFile\"}"."\n";
	$seqin = new Bio::SeqIO('-file'=>$self->{"readsFile"});
    }
    while (my $seq = $seqin->next_seq) {
	if(exists $hits{$seq->id}){
	    foreach my $markerHit(keys %{$hits{$seq->id}}){
		#checking if a 6frame translation was done and the suffix was appended to the description and not the sequence ID
		my $newID = $seq->id;
		my $current_suff="";
		my $current_seq="";
		if(exists $duplicates{$current_seq}{$markerHit}{$current_suff}){
		    warn "Skipping ".$seq->id."\t".$current_suff."\n";
		    next;
		}

		#pre-trimming for the query + FLANKING_LENGTH residues before and after (for very long queries)
		my $start = $hitsStart{$seq->id}{$markerHit};
		my $end = $hitsEnd{$seq->id}{$markerHit};
		($start,$end) = ($end,$start) if($start > $end); # swap if start bigger than end
		$start -= FLANKING_LENGTH;
		$end += FLANKING_LENGTH;
		$start=abs($start) % 3 + 1 if($start < 0);
		my $seqLength = length($seq->seq);
		$end= $end-ceil(($end - $seqLength)/3)*3 if($end >= $seqLength);

		my $newSeq = substr($seq->seq,$start,$end-$start);
		#if the $newID exists in the %frames hash, then it needs to be translated to the correct frame
		if($reverseTranslate){
		    # compute the frame as modulo 3 of start site, reverse strand if end < start
		    my $frame = $hitsStart{$seq->id}{$markerHit} % 3 + 1;
		    $frame *= -1 if( $hitsStart{$seq->id}{$markerHit} > $hitsEnd{$seq->id}{$markerHit});
		    my $seqlen = abs($hitsStart{$seq->id}{$markerHit} - $hitsEnd{$seq->id}{$markerHit})+1;
		    if($seqlen % 3 == 0){
	                    $newSeq = translateFrame($newID,$seq->seq,$start,$end,$frame,$markerHit,$reverseTranslate);
		    }else{
			warn "Error, alignment length not multiple of 3!  FIXME: need to pull frameshift from full blastx\n";
		    }
                }
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
	open(fileOUT,">".$self->{"blastDir"}."/$marker.candidate")or die " Couldn't open ".$self->{"blastDir"}."/$marker.candidate for writing\n";
	print fileOUT $markerHits{$marker};
	close(fileOUT);
	if($reverseTranslate){
	    open(fileOUT,">".$self->{"blastDir"}."/$marker.candidate.ffn")or die " Couldn't open ".$self->{"blastDir"}."/$marker.candidate.ffn for writing\n";
	    print fileOUT $markerNuc{$marker};
	    close(fileOUT);
	}
	
    }

    return $self;
}

=head2 prepAndClean

=item *

Checks if the directories needed for the blast run and parsing exist
Removes previous blast runs data if they are still in the directories
Generates the blastable database using the marker representatives

=back

=cut

sub prepAndClean {
    my $self = shift;
    debug "prepclean MARKERS @markers\nTESTING\n ";
    `mkdir $self->{"tempDir"}` unless (-e $self->{"tempDir"});
    #create a directory for the Reads file being processed.
    `mkdir $self->{"fileDir"}` unless (-e $self->{"fileDir"});
    `mkdir $self->{"blastDir"}` unless (-e $self->{"blastDir"});
    if($custom ne ""){
	#remove rep.faa if it already exists (starts clean is a previous job was stopped or crashed or included different markers)
	if(-e $self->{"blastDir"}."/rep.faa"){ `rm $self->{"blastDir"}/rep.faa`;}
	#also makes 1 large file with all the marker sequences
	open( REPFILE, ">".$self->{"blastDir"}."/rep.faa" );
	foreach my $marker (@markers){
	    #if a marker candidate file exists remove it, it is from a previous run and could not be related
	    if(-e $self->{"blastDir"}."/$marker.candidate"){
		`rm $self->{"blastDir"}/$marker.candidate`;
	    }
	    #initiate the hash table for all markers incase 1 marker doesn't have a single hit, it'll still be in the results 
	    #and will yield an empty candidate file
	    $markerHits{$marker}="";
	    #append the rep sequences for all the markers included in the study to the rep.faa file
	    my $fastaName = Amphora2::Utilities::getFastaMarkerFile($self,$marker);
	    open( FASTA, "$Amphora2::Utilities::marker_dir/$fastaName" );
		while( my $fline = <FASTA> ){
			if($fline =~ /^>/){
				my $mcat = $marker."__";
				$fline =~ s/^>/>$mcat/g;
			}
			print REPFILE $fline;
		}
#	    `cat $Amphora2::Utilities::marker_dir/$fastaName.faa >> $self->{"blastDir"}/rep.faa`
	}
	close REPFILE;
	#make a DB for rapSearch
	if(!-e $self->{"blastDir"}."/rep.des" ||  !-e $self->{"blastDir"}."/rep.fbn" || !-e $self->{"blastDir"}."/rep.inf" || !-e $self->{"blastDir"}."/rep.swt"){
	    `cd $self->{"blastDir"} ; $Amphora2::Utilities::preRapSearch -d $self->{"blastDir"}/rep.faa -n rep`;
        }
	#make one for BLAST
	if(!-e $self->{"blastDir"}."/rep.faa.psq" ||  !-e $self->{"blastDir"}."/rep.faa.pin" || !-e $self->{"blastDir"}."/rep.faa.phr"){
	    `formatdb -i $self->{"blastDir"}/rep.faa -o F -p T -t RepDB`;
        }
    }else{
	#when using the default marker package
	my $dbDir = "$Amphora2::Utilities::marker_dir/representatives";
	debug "Using the standard marker package\n";
	debug "Using $dbDir as default directory\n";
	#removing all the previous representative work from other runs in case we switch from Updated to stock markers
	if(-e "$dbDir/rep.faa"){
	    `rm $dbDir/rep.*`;
	}
	if(!-e "$dbDir/rep.faa"){
		open( REPFILE, ">$dbDir/rep.faa" );
	    foreach my $marker (@markers){
		$markerHits{$marker}="";
		    my $fastaName = Amphora2::Utilities::getFastaMarkerFile($self,$marker);
		    open( FASTA, "$Amphora2::Utilities::marker_dir/$fastaName" );
			while( my $fline = <FASTA> ){
				if($fline =~ /^>/){
					my $mcat = $marker."__";
					$fline =~ s/^>/>$mcat/g;
				}
				print REPFILE $fline;
			}
		}
		close REPFILE;
	}
	if(!-e "$dbDir/rep.des" ||  !-e "$dbDir/rep.fbn" || !-e "$dbDir/rep.inf" || !-e "$dbDir/rep.swt"){
		`cd $dbDir ; $Amphora2::Utilities::preRapSearch -d rep.faa -n rep`;
	}
	# now make a DB for BLAST containing just representative sequences
	if(!-e "$dbDir/$blastdb_name"){
	    foreach my $marker (@markers){
		`cat $dbDir/$marker.rep >> $dbDir/$blastdb_name`;
	    }
	}
    `cp $dbDir/$blastdb_name $self->{"blastDir"}`;
    `$Amphora2::Utilities::formatdb -i $self->{"blastDir"}/$blastdb_name -o F -p T -t RepDB`;
    }
    return $self;
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

    perldoc Amphora2::blast


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

1; # End of Amphora2::blast.pm
