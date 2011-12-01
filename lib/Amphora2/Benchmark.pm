package Amphora2::Benchmark;

use warnings;
use strict;
use Carp;
use Bio::Phylo;
use Amphora2::Summarize qw(:all);

=head1 SUBROUTINES/METHODS

=Head2 benchmark_illumina

WARNING: The input file must have the correct format to run the benchmark
Reads the input file looking for taxonomic origins in the read headers and 
Reads Summary files and generates accuracy ratings for every taxonmic level.

=cut


my %parent=();
my %nameidmap=();  # Key is an id - Value is a name (Use if you have an ID and want a name)
my %idnamemap=();  # Key is a name - Value is an id (Use if you have a name and want an ID)
my %sourceIDs=();
my %readSource=();
my %correctReadPlacement=();
my %refTaxa = ();
sub runBenchmark{
    my $self = shift;
    my $output_path = shift;
    (%nameidmap, %idnamemap)= Amphora2::Summarize::readNcbiTaxonNameMap();
    %parent = Amphora2::Summarize::readNcbiTaxonomyStructure();
    %refTaxa = getInputTaxa($self->{"readsFile"});
    print "Number of reads counted = ".scalar(keys(%readSource))."\n";
    
    readSeqSummary($self, $output_path, \%readSource);

}

=head2 readSeqSummary
Takes a Directory name for input
Reads the sequence_taxa.txt file and compares the read placements to their true source
prints the percentage of all reads that have the correct taxanomic ID
prints the percentage of all PLACED reads that have the correct taxonmic ID
=cut

sub readSeqSummary{
    my $self = shift;
    my $output_path = shift;
    my $readSource = shift;
    my $targetDir = $self->{"fileDir"};
    open(fileIN,$targetDir."/sequence_taxa.txt");
    my %topReadScore = (); 
    my %allPlacedScore = ();
    #reading and storing information from the sequence_taxa.txt file
    while(<fileIN>){
	if($_ =~ m/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/){
	    my $read = $1;
	    my $taxPlacement = $4;
	    my $probability = $5;
	    my @taxPlacementID = Amphora2::Summarize::getTaxonInfo($2);
#	    print "TaxPlacement : $2\t $taxPlacementID[0]\t\n";

	    my @readAncestor = getAncestorArray($taxPlacementID[2]);
#	    print $read." $taxPlacementID[1]\t$taxPlacementID[2]:\t@readAncestor\n";
	    my $rank = $taxPlacementID[1];
	    #keep only the top hits for all ranks for each Read
	    if(exists $topReadScore{$read}){
		if($topReadScore{$read}->[0] < $probability ){
		    if( $topReadScore{$read}->[2] < scalar(@readAncestor)){
			my @array = ($probability,$taxPlacementID[2],scalar(@readAncestor));
			$topReadScore{$read} = \@array;
		    }
		}else{
		    # Do nothing
		}
	    }else{
		my @array =($probability,$taxPlacementID[2],scalar(@readAncestor));
		$topReadScore{$read} = \@array;
	    }
	    my @array =($probability,$taxPlacementID[2],scalar(@readAncestor));
	    $allPlacedScore{$read}{$taxPlacement} = \@array;
#	    exit;
	}else{
#	    debug "Warning : Line format not recognized .... Skipping line from sequence_taxa.txt\n";
	}
    }
    close(fileIN);

    #comparing the sequence_taxa information with the Source taxons
#    my %overallScore;
    my %matchTop=();
    foreach my $readID (keys %topReadScore){
	#look at each taxonomic level for each Read
	my @ancArrayRead = getAncestorArray($topReadScore{$readID}->[1]);
	pop(@ancArrayRead);
	push(@ancArrayRead,$topReadScore{$readID}->[1]);
	foreach my $id (@ancArrayRead){
	    if(exists $sourceIDs{$id}){
		my @currTaxon = Amphora2::Summarize::getTaxonInfo($id);
		my $currRank = $currTaxon[1];
		if(exists $matchTop{$currRank}){
		    $matchTop{$currRank}++;
		}else{
		    $matchTop{$currRank}=1;
		}
	    }
	}
    }
    my $readNumber = scalar(keys(%topReadScore));

    my $allReadNumber = 0;
    my $totalProb=0;
    my %rankTotalProb = ();
    my %matchAll=();
    foreach my $readID (keys %allPlacedScore){
	foreach my $tax (keys %{$allPlacedScore{$readID}}){
	    $allReadNumber++;
	    my @ancArrayRead = getAncestorArray($allPlacedScore{$readID}{$tax}->[1]);
	    pop(@ancArrayRead);
	    push(@ancArrayRead,$allPlacedScore{$readID}{$tax}->[1]);
	    foreach my $id (@ancArrayRead){
		my @currTaxon = Amphora2::Summarize::getTaxonInfo($id);
		my $currRank = $currTaxon[1];
		if(exists $sourceIDs{$id}){
		    if(exists $matchAll{$currRank}){
			$matchAll{$currRank} += $allPlacedScore{$readID}{$tax}->[0];
		    }else{
			$matchAll{$currRank}=$allPlacedScore{$readID}{$tax}->[0];
		    }
		}
		if(exists $rankTotalProb{$currRank}){
		    $rankTotalProb{$currRank} += $allPlacedScore{$readID}{$tax}->[0];
		}else{
		    $rankTotalProb{$currRank} = $allPlacedScore{$readID}{$tax}->[0];
		}
		$totalProb +=$allPlacedScore{$readID}{$tax}->[0];
	    }
        }
    }
	chdir($output_path) if defined($output_path);
	report_text($self, "accuracy.txt", \%matchTop, \%matchAll, $readNumber, $allReadNumber, $totalProb, \%rankTotalProb);
	report_csv($self, \%matchTop, \%matchAll, $readNumber, $allReadNumber, $totalProb, \%rankTotalProb);
}

#    foreach my $m (keys %matchTop){
#	next if $m eq "no rank";
	#print "Match : $m\t".$match{$m}/$readNumber."\n";
#	print "Top placed Matches : $m\t".$matchTop{$m}."\n";
#}

sub report_flot_json{
	my $mtref = shift;
	my $maref = shift;
	my $readNumber = shift;
	my $allReadNumber = shift;
	my $totalProb = shift;
	my $rtpref = shift;
	my %matchTop = %$mtref;
	my %matchAll = %$maref;
	my %rankTotalProb = %$rtpref;
	
#	open(JSON, ">amphora2_accuracy.json");
}

# get the date in YYYYMMDD format
sub get_date_YYYYMMDD {
	my @timerval = localtime();
	my $datestr = (1900+$timerval[5]);
	$datestr .= 0 if $timerval[4] < 9; 
	$datestr .= ($timerval[4]+1);
	$datestr .= 0 if $timerval[3] < 9; 
	$datestr .= $timerval[3];
	return $datestr;
}

sub reportTiming{
	my $self = shift;
	my $data = shift;
	my $output_path = shift;
	my $timing_file = $output_path."/timing.csv";
	unless(-f $timing_file){
		open(TIMING, ">$timing_file");
		print TIMING "Date,".join(",",keys(%$data))."\n";
		close TIMING;
	}
	open(TIMING, ">>$timing_file");
	print TIMING get_date_YYYYMMDD;
	foreach my $time(keys(%$data)){
		print TIMING ",".$data->{$time};
	}
	print TIMING "\n";
}

sub report_csv{
	my $self = shift;
	my $mtref = shift;
	my $maref = shift;
	my $readNumber = shift;
	my $allReadNumber = shift;
	my $totalProb = shift;
	my $rtpref = shift;
	my %matchTop = %$mtref;
	my %matchAll = %$maref;
	my %rankTotalProb = %$rtpref;

	my $tophitfile = $self->{"readsFile"}."tophits.csv";
	unless(-f $tophitfile){
		open(TOPHITS, ">$tophitfile");
		print TOPHITS "Date,Superkingdom,Phylum,Subphylum,Class,Order,Family,Genus,Species,Subspecies,No Rank\n";
		close TOPHITS;
	}

	my $date = get_date_YYYYMMDD();

	# append an entry to the tophits file
	open(TOPHITS, ">>$tophitfile");
	print TOPHITS $date;
	print TOPHITS ",".100*$matchTop{"superkingdom"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"phylum"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"subphylum"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"class"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"order"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"family"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"genus"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"species"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"subspecies"}/$readNumber;
	print TOPHITS ",".100*$matchTop{"no rank"}/$readNumber;
	print TOPHITS "\n";
	
}

sub report_text{
    my $outputfile = shift;
    my $mtref = shift;
    my $maref = shift;
    my $readNumber = shift;
    my $allReadNumber = shift;
    my $totalProb = shift;
    my $rtpref = shift;
    my %matchTop = %$mtref;
    my %matchAll = %$maref;
    my %rankTotalProb = %$rtpref;

    print "\n";
    print "Top placed Matches : Superkingdom\t".100*$matchTop{"superkingdom"}/$readNumber."\n";
    print "Top placed Matches : Phylum\t".100*$matchTop{"phylum"}/$readNumber."\n";
    print "Top placed Matches : Subphylum\t".100*$matchTop{"subphylum"}/$readNumber."\n";
    print "Top placed Matches : Class\t".100*$matchTop{"class"}/$readNumber."\n";
    print "Top placed Matches : Order\t".100*$matchTop{"order"}/$readNumber."\n";
    print "Top placed Matches : Family\t".100*$matchTop{"family"}/$readNumber."\n";
    print "Top placed Matches : Genus\t".100*$matchTop{"genus"}/$readNumber."\n";
    print "Top placed Matches : Species\t".100*$matchTop{"species"}/$readNumber."\n";
    print "Top placed Matches : Subspecies\t".100*$matchTop{"subspecies"}/$readNumber."\n";
    print "Top placed Matches : No rank\t".100*$matchTop{"no rank"}/$readNumber."\n";
    print "\n";
    print "Total placed Reads : $readNumber\n";

    print "\n";
    print "All placements Matches : Superkingdom\t".100*$matchAll{"superkingdom"}/$allReadNumber."\n";
    print "All placements Matches : Phylum\t".100*$matchAll{"phylum"}/$allReadNumber."\n";
    print "All placements Matches : Subphylum\t".100*$matchAll{"subphylum"}/$allReadNumber."\n";
    print "All placements Matches : Class\t".100*$matchAll{"class"}/$allReadNumber."\n";
    print "All placements Matches : Order\t".100*$matchAll{"order"}/$allReadNumber."\n";
    print "All placements Matches : Family\t".100*$matchAll{"family"}/$allReadNumber."\n";
    print "All placements Matches : Genus\t".100*$matchAll{"genus"}/$allReadNumber."\n";
    print "All placements Matches : Species\t".100*$matchAll{"species"}/$allReadNumber."\n";
    print "All placements Matches : Subspecies\t".100*$matchAll{"subspecies"}/$allReadNumber."\n";
    print "All placements Matches : No rank\t".100*$matchAll{"no rank"}/$allReadNumber."\n";
    print "\n";
    print "All placements Matches : Superkingdom\t".100*$matchAll{"superkingdom"}/$totalProb."\n";
    print "All placements Matches : Phylum\t".100*$matchAll{"phylum"}/$totalProb."\n";
    print "All placements Matches : Subphylum\t".100*$matchAll{"subphylum"}/$totalProb."\n";
    print "All placements Matches : Class\t".100*$matchAll{"class"}/$totalProb."\n";
    print "All placements Matches : Order\t".100*$matchAll{"order"}/$totalProb."\n";
    print "All placements Matches : Family\t".100*$matchAll{"family"}/$totalProb."\n";
    print "All placements Matches : Genus\t".100*$matchAll{"genus"}/$totalProb."\n";
    print "All placements Matches : Species\t".100*$matchAll{"species"}/$totalProb."\n";
    print "All placements Matches : Subspecies\t".100*$matchAll{"subspecies"}/$totalProb."\n";
    print "All placements Matches : No rank\t".100*$matchAll{"no rank"}/$totalProb."\n";
    print "\n";
    print "Rank specific Percentages\n";
    print "All placements Matches : Superkingdom\t".100*$matchAll{"superkingdom"}/$rankTotalProb{"superkingdom"}."\n";
    print "All placements Matches : Phylum\t".100*$matchAll{"phylum"}/$rankTotalProb{"phylum"}."\n";
    print "All placements Matches : Subphylum\t".100*$matchAll{"subphylum"}/$rankTotalProb{"subphylum"}."\n";
    print "All placements Matches : Class\t".100*$matchAll{"class"}/$rankTotalProb{"class"}."\n";
    print "All placements Matches : Order\t".100*$matchAll{"order"}/$rankTotalProb{"order"}."\n";
    print "All placements Matches : Family\t".100*$matchAll{"family"}/$rankTotalProb{"family"}."\n";
    print "All placements Matches : Genus\t".100*$matchAll{"genus"}/$rankTotalProb{"genus"}."\n";
    print "All placements Matches : Species\t".100*$matchAll{"species"}/$rankTotalProb{"species"}."\n";
    print "All placements Matches : Subspecies\t".100*$matchAll{"subspecies"}/$rankTotalProb{"subspecies"}."\n";
    print "All placements Matches : No rank\t".100*$matchAll{"no rank"}/$rankTotalProb{"no rank"}."\n";
    print "\n";





    print "Total placements : $allReadNumber\n";
}

=head2 getInputTaxa
Reads a fasta file extracting the source field for each read.
Compiles statistics on the input read abundance
Returns a hash of taxonomic ancestry for Source organisms

TODO : Determine which reads came from the marker gene regions from the source genomes

=cut

sub getInputTaxa{
    my $fileName = shift;
    my %sourceTaxa=();
    my %sourceReadCounts=();
    open(fileIN,$fileName) or carp("Couldn't open ".$fileName."\n");
    while(<fileIN>){
	next unless $_=~m/^>/;
	$_ =~ m/^>(\S+).*SOURCE_\d+="(.*)"/;
	#push(@sourceTaxa,$1);
	$readSource{$1}=$2;
	my @ancestors = getAncestorArray($2);
	foreach my $id (@ancestors){
	    $sourceIDs{$id}=1;
	}
	$sourceIDs{$2}=1;
	if(exists $sourceReadCounts{$2}){
	    $sourceReadCounts{$2}++;
	}else{
	    $sourceReadCounts{$2}=0;
	    $sourceTaxa{$2}=getAncestorArray($2);
	}
    }
    close(fileIN);

    foreach my $source(keys %sourceReadCounts){
	print $source."\t".$nameidmap{$source}."\t".$sourceReadCounts{$source}."\n";
    }

    return %sourceTaxa;
}

=head2 getAncestorArray
Takes an input taxID and returns the ancestors in an array going up the tree.
index 0 being the first ancestor.
last index being the root of the tree.
=cut

sub getAncestorArray{
    my $taxID = shift;
    my $curID = $taxID;
    my @ancestor=();
    while($curID !=1){
	$curID = ${$parent{$curID}}[0];
	push(@ancestor,$curID);
    }
    return @ancestor;
}

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Amphora2::Benchmark


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

1; # End of Amphora2::Benchmark.pm
