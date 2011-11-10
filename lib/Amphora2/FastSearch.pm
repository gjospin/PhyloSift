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

    # ensure databases and sequences are prepared for search
    debug "before rapPrepandclean\n";
    prepAndClean($self);
    if($self->{"readsFile_2"} ne ""){
	debug "before fastqtoFASTA\n";
	fastqToFasta($self);
    }

    # search reads/contigs against marker database
    my $resultsfile;
    if(!defined($self->{"dna"}) || $self->{"dna"}==0){
	$resultsfile = executeBlast($self, $self->{"readsFile"});
    }elsif($searchtype eq "blast"){
	$resultsfile = blastXoof_table($self, $self->{"readsFile"});
	$reverseTranslate=1;
    }else{
        $resultsfile = executeRap($self);
	build_lookup_table($self);
    }

    # parse the hits to marker genes
    my $hitsref;
    if(defined $self->{"coverage"} && (!defined $self->{"isolate"} || $self->{"isolate"} != 1)){
	$hitsref = get_hits_contigs($self,$resultsfile, $searchtype);
    }else{
	$hitsref = get_hits($self,$resultsfile, $searchtype);
    }

    # write out sequence regions hitting marker genes to candidate files
    writeCandidates($self,$hitsref);
    
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

=head2 get_hits_contigs

parse the blast file

=cut


sub get_hits_contigs{
	my $self = shift;
	my $hitfilename=shift;

	# key is a contig name
	# value is an array of arrays, each one has [marker,bit_score,left-end,right-end]
	my %contig_hits;
	my %contig_top_bitscore;
	my $max_hit_overlap = 10;

	open(blastIN,$hitfilename)or carp("Couldn't open ".$hitfilename."\n");
	while(<blastIN>){
		# read a blast line
		next if($_ =~ /^#/);
		chomp($_);
		my ($query, $subject, $two, $three, $four, $five, $query_start,
			$query_end, $eight, $nine, $ten, $bitScore) = split(/\t/,$_);	
		# get the marker name
		my @marker=split(/\_/,$subject);
		my $markerName = $marker[$#marker];
		# running on long reads or an assembly
		# allow each region of a sequence to have a top hit
		# do not allow overlap
		if(defined($contig_top_bitscore{$query}{$markerName}) && $contig_top_bitscore{$query}{$markerName} - $bestHitsBitScoreRange < $bitScore){
			my $i=0;
			for(; $i<@{$contig_hits{$query}}; $i++){
				my $prevhitref=$contig_hits{$query}->[$i];
				my @prevhit = @$prevhitref;
				# is there enough overlap to consider these the same?
				# if so, take the new one if it has higher bitscore
				if(	$prevhit[2] < $prevhit[3] && $query_start < $query_end &&
					$prevhit[2] < $query_end - $max_hit_overlap && 
					$query_start + $max_hit_overlap < $prevhit[3]){
					print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [$markerName, $bitScore, $query_start, $query_end] if ($bitScore > $prevhit[1]);
					last;
				}
				# now check the same for reverse-strand hits
				if(	$prevhit[2] > $prevhit[3] && $query_start > $query_end &&
					$prevhit[3] < $query_start - $max_hit_overlap && 
					$query_end + $max_hit_overlap < $prevhit[2]){
					print STDERR "Found overlap $query and $markerName, $query_start:$query_end\n";
					$contig_hits{$query}->[$i] = [$markerName, $bitScore, $query_start, $query_end] if ($bitScore > $prevhit[1]);
					last;
				}
			}
			if($i==@{$contig_hits{$query}}){
				# no overlap was found, include this hit
				my @hitdata = [$markerName, $bitScore, $query_start, $query_end];
				push(@{$contig_hits{$query}}, @hitdata);
			}
		}elsif(!defined($contig_top_bitscore{$query}{$markerName})){
			my @hitdata = [$markerName, $bitScore, $query_start, $query_end];
			push(@{$contig_hits{$query}}, @hitdata);
			$contig_top_bitscore{$query}{$markerName} = $bitScore;
		}
	}
	return \%contig_hits;
}


=head2 get_hits

parse the blast file

=cut

sub get_hits{
    my $self = shift;
    my $hitfilename=shift;
    my $searchtype = shift;

    my %markerTopScores;
    my %topScore=();
    my %contig_hits;

    open(blastIN,$hitfilename)or carp("Couldn't open ".$hitfilename."\n");
    while(<blastIN>){
	chomp($_);
	next if($_ =~ /^#/);
	my ($query, $subject, $two, $three, $four, $five, $query_start,
		$query_end, $eight, $nine, $ten, $bitScore) = split(/\t/,$_);	
	my $markerName = getMarkerName($subject, $searchtype);

#	if($searchtype ne "blast" && $self->{"dna"}){
		# RAPsearch seems to have an off-by-one on its DNA coordinates
#		$query_start -= 2;
#		$query_end -= 2;
#	}
	#parse once to get the top score for each marker (if isolate is ON, parse again to check the bitscore ranges)
	if($isolateMode==1){
	    # running on a genome assembly, allow only 1 hit per marker (TOP hit)
	    if( !defined($markerTopScores{$markerName}) || $markerTopScores{$markerName} < $bitScore ){
		$markerTopScores{$markerName} = $bitScore;
	    }
	}else{
	    # running on short reads, just do one marker per read
		$topScore{$query}=0 unless exists $topScore{$query};
		#only keep the top hit
		if($topScore{$query} <= $bitScore){
		    $contig_hits{$query} = [[$markerName, $bitScore, $query_start, $query_end]];
		    $topScore{$query}=$bitScore;
		}#else do nothing
	}
    }
    close(blastIN);
    if($isolateMode==1){
	# reading the output a second to check the bitscore ranges from the top score
	open(blastIN,$hitfilename)or die "Couldn't open $hitfilename\n";
	# running on a genome assembly, allow more than one marker per sequence
	# require all hits to the marker to have bit score within some range of the top hit
	while(<blastIN>){
		chomp($_);
		next if($_ =~ /^#/);
		my ($query, $subject, $two, $three, $four, $five, $query_start,
			$query_end, $eight, $nine, $ten, $bitScore) = split(/\t/,$_);	
		my $markerName = getMarkerName($subject, $searchtype);
		if($markerTopScores{$markerName} < $bitScore + $bestHitsBitScoreRange){
			my @hitdata = [$markerName, $bitScore, $query_start, $query_end];
			push(@{$contig_hits{$query}}, @hitdata);
		}
	}
	close(blastIN);
    }
    return \%contig_hits;
}


=head2 getMarkerName

Extracts a marker gene name from a blast or rapsearch subject sequence name

=cut

sub getMarkerName{
	my $subject = shift;
	my $searchtype = shift;
	my $markerName = "";
	if($searchtype eq "blast"){
		my @marker=split(/\_/,$subject);
		$markerName = $marker[$#marker];
	}else{
		my @marker=split(/\_\_/,$subject);
		$markerName = $marker[0];
#		$markerName = $marker_lookup{$subject};
	}
	return $markerName;
}

=head2 writeCandidates

write out results

=cut

sub writeCandidates{
	my $self = shift;
	my $contigHitsRef = shift;
	my %contig_hits = %$contigHitsRef;
	debug "ReadsFile:  $self->{\"readsFile\"}"."\n";
	my $seqin = new Bio::SeqIO('-file'=>$self->{"readsFile"});
	while (my $seq = $seqin->next_seq) {
		# skip this one if there are no hits
		next unless( exists $contig_hits{$seq->id} );
		for( my $i=0; $i<@{$contig_hits{$seq->id}}; $i++){
			my $curhitref=$contig_hits{$seq->id}->[$i];
			my @curhit=@$curhitref;
			my $markerHit = $curhit[0];
			my $start = $curhit[2];
			my $end = $curhit[3];
			($start,$end) = ($end,$start) if($start > $end); # swap if start bigger than end
			$start -= FLANKING_LENGTH;
			$end += FLANKING_LENGTH;
			# ensure flanking region is a multiple of 3 to avoid breaking frame in DNA
			$start=abs($start) % 3 + 1 if($start < 0);
			my $seqLength = length($seq->seq);
			$end= $end-ceil(($end - $seqLength)/3)*3 if($end >= $seqLength);

			my $newSeq = substr($seq->seq,$start,$end-$start);
			#if we're working from DNA then need to translate to protein
			if($self->{"dna"}){
				# compute the frame as modulo 3 of start site, reverse strand if end < start
				my $frame = $curhit[2] % 3 + 1;
				$frame *= -1 if( $curhit[2] > $curhit[3]);
				my $seqlen = abs($curhit[2] - $curhit[3])+1;
				if($seqlen % 3 == 0){
					$newSeq = translateFrame($seq->id,$seq->seq,$start,$end,$frame,$markerHit,$self->{"dna"});
					$newSeq =~ s/\*/X/g;	# bioperl uses * for stop codons but we want to give X to hmmer later
				}else{
					warn "Error, alignment length not multiple of 3!  FIXME: need to pull frameshift from full blastx\n";
				}
			}
			$markerHits{$markerHit} = "" unless defined($markerHits{$markerHit});
			$markerHits{$markerHit} .= ">".$seq->id."\n".$newSeq."\n";
#			$markerHits{$markerHit} .= ">".$seq->id.":$start-$end\n".$newSeq."\n";
		}		
	}

	#write the read+ref_seqs for each markers in the list
	foreach my $marker (keys %markerHits){
		#writing the hits to the candidate file
		open(fileOUT,">".$self->{"blastDir"}."/$marker.candidate")or die " Couldn't open ".$self->{"blastDir"}."/$marker.candidate for writing\n";
		print fileOUT $markerHits{$marker};
		close(fileOUT);
		if($self->{"dna"}){
			open(fileOUT,">".$self->{"blastDir"}."/$marker.candidate.ffn")or die " Couldn't open ".$self->{"blastDir"}."/$marker.candidate.ffn for writing\n";
			print fileOUT $markerNuc{$marker};
			close(fileOUT);
		}
	}
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
