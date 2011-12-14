package Amphora2::Summarize;

use warnings;
use strict;
use FindBin;
use Amphora2::Amphora2;
use Carp;
use Bio::Phylo;
use Bio::Phylo::Forest::Tree;
if($^O=~/arwin/){
	use lib "$FindBin::Bin/osx/darwin-thread-multi-2level/";
}

=head1 NAME

Amphora2::Summarize - Summarize placed reads using the NCBI taxonomy

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Reconciles gene-tree specific read placements with the NCBI taxonomy

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 summarize

=cut

my %nameidmap;
my %idnamemap;

# read the NCBI taxon names
# stash them in hashes called nameidmap and idnamemap to go back & forth from tax ids to names
sub readNcbiTaxonNameMap {
    my $ncbidir = $Amphora2::Utilities::ncbi_dir;
    open( TAXIDS, "$ncbidir/names.dmp" );
    while( my $line = <TAXIDS> ){
	chomp $line;
	if(($line =~ /scientific name/) || ($line =~ /synonym/) || ($line =~ /misspelling/)){
	    my @vals = split( /\s+\|\s+/, $line );
	    $nameidmap{homogenizeNameAlaDongying($vals[1])}=$vals[0];
	    $idnamemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]) if($line =~ /scientific name/);
	}
    }
    return (%nameidmap,%idnamemap);
}

# now read the NCBI taxonomy structure
# puts the results in a hash called "parent"
my %parent;
sub readNcbiTaxonomyStructure {
    my $ncbidir = $Amphora2::Utilities::ncbi_dir;
    open( TAXSTRUCTURE, "$ncbidir/nodes.dmp" );
    while( my $line = <TAXSTRUCTURE> ){
	chomp $line;
	my @vals = split( /\s+\|\s+/, $line );
	$parent{$vals[0]} = [$vals[1],$vals[2]];
    }
    return %parent;
}


sub makeNcbiTreeFromUpdate {
	my $self = shift;
	my $results_dir = shift;
	my $markerdir = shift;
	readNcbiTaxonNameMap();
	readNcbiTaxonomyStructure();
	my @orgnames = `ls -1 $results_dir | grep fasta`; 
	my @taxonids;
	open( MARKERTAXONMAP, ">$markerdir/marker_taxon_map.updated.txt" );
	foreach my $org(@orgnames){
		$org =~ /\.(\d+)\.fasta/;
		if(!defined($1)){
			print MARKERTAXONMAP treeName($org)."\t".treeName($org)."\n";
			next;
		}
		print STDERR "Bad taxon ID $1 for $org" unless defined($parent{$1});
		next unless defined($parent{$1});
		push(@taxonids, $1);
		chomp($org);
		print MARKERTAXONMAP "$1\t$1\n";
	}
	close MARKERTAXONMAP;
	my %tidnodes;
	my $phylotree = Bio::Phylo::Forest::Tree->new();
	foreach my $tid(@taxonids){
		next if ($tid eq "");
		my $child;        
		while( $tid != 1 ){
		    # check if we've already seen this one
		    last if(defined($tidnodes{$tid}));
		    # create a new node & add to tree
		    my $parentid = $parent{$tid}->[0];
		    my $newnode;
		    $newnode = Bio::Phylo::Forest::Node->new( -parent=>$tidnodes{$parentid}, -name=>$tid) if defined($tidnodes{$parentid});
		    $newnode = Bio::Phylo::Forest::Node->new( -name=>$tid) if !defined($tidnodes{$parentid});
		    $tidnodes{$tid} = $newnode;
		    $newnode->set_child($child) if(defined($child));
		    $phylotree->insert($newnode);
		    # continue traversal toward root
		    $tid = $parentid;
		    $child = $newnode;
		}
	}
	open( TREEOUT, ">ncbi_tree.updated.tre" );
	print TREEOUT $phylotree->to_newick("-nodelabels"=>1);
	close TREEOUT;
}

=head2 makeNcbiTree
Reads all the marker gene trees, finds their corresponding taxa in the NCBI taxonomy, and
constructs a newick format tree representation of the NCBI taxonomy containing only the
organisms present in the marker gene trees.  
=cut
sub makeNcbiTree {
    my $self = shift;
    readNcbiTaxonNameMap();
    readNcbiTaxonomyStructure();
    # now read the list of organisms we have in our DB
    # construct a phylo tree with the NCBI topology containing
    # just the organisms in our database
    my $markerdir = $Amphora2::Utilities::marker_dir;
    my %namemap = Amphora2::Utilities::readNameTable($markerdir);
    my $phylotree = Bio::Phylo::Forest::Tree->new();
    open( MARKERTAXONMAP, ">$markerdir/marker_taxon_map.txt" );
    my %tidnodes;
    foreach my $key( keys(%namemap) ){
	$namemap{$key}=homogenizeNameAlaDongying($namemap{$key});
	my ($tid,$name) = dongyingFindNameInTaxaDb($namemap{$key});
	if($tid eq "ERROR"){
	    print STDERR "Error! Could not find $namemap{$key} in name map\n" if length($key) > 12;
	    next;
	}
	# add it to the mapping file
	my $treename = treeName($idnamemap{$tid});
        print MARKERTAXONMAP "$key\t$treename\n";
	#got the taxon id, now walk to root adding tree nodes as necessary
	next unless(defined($tid));
	my $child;        
	while( $tid != 1 ){
	    # check if we've already seen this one
            last if(defined($tidnodes{$tid}));
            # create a new node & add to tree
            my $nodename = treeName($idnamemap{$tid});
	    my $parentid = $parent{$tid}->[0];
            my $newnode;
            $newnode = Bio::Phylo::Forest::Node->new( -parent=>$tidnodes{$parentid}, -name=>$nodename) if defined($tidnodes{$parentid});
            $newnode = Bio::Phylo::Forest::Node->new( -name=>$nodename) if !defined($tidnodes{$parentid});
            $tidnodes{$tid} = $newnode;
            $newnode->set_child($child) if(defined($child));
            $phylotree->insert($newnode);
            # continue traversal toward root
            $tid = $parentid;
            $child = $newnode;
	}
    }
    close MARKERTAXONMAP;
    open( TREEOUT, ">ncbi_tree.tre" );
    print TREEOUT $phylotree->to_newick("-nodelabels"=>1);
    close TREEOUT;
}

=head2 readCoverage
Reads a coverage file
=cut
sub readCoverage {
	my $coverageFile = shift;
	my %coverage;
	open( COVERAGE, $coverageFile ) || return %coverage;
	while( my $line = <COVERAGE> ){
		chomp $line;
		my @data = split( /\t/, $line );
		$data[0] =~ s/[\(\)\[\]\+\=\<\>\?]//g;
		$data[0] =~ s/[\-\s]/_/g;
		$coverage{$data[0]} = $data[1];
	}
	return \%coverage;
}

=head2 summarize
Reads the .place files containing Pplacer read placements and maps them onto the
NCBI taxonomy
=cut
sub summarize {
    my $self = shift;
    my $markRef = shift; # list of the markers we're using
    readNcbiTaxonNameMap();
    readNcbiTaxonomyStructure();
    my $markerdir = $Amphora2::Utilities::marker_dir;
    my %namemap = Amphora2::Utilities::readNameTable($markerdir);
    foreach my $key( keys(%namemap) ){
	$namemap{$key}=homogenizeNameAlaDongying($namemap{$key});
    }
    # keep a hash counting up all the read placements
    my %ncbireads;

    # try to read a contig coverage file if it exists
    my %coverage;

    if(defined $self->{"coverage"} ){
	my $covref = readCoverage($self->{"coverage"});
	%coverage = %$covref;
    }

    # read all of the .place files for markers
    # map them onto the ncbi taxonomy
    # this is a hash structured as {sequenceID}{taxonID}=probabilitySum
    my %placements;
    foreach my $marker(@{$markRef}){
	# don't bother with this one if there's no read placements
	my $placeFile = $self->{"treeDir"}."/".Amphora2::Utilities::getReadPlacementFile($marker);
	next unless( -e $placeFile );

        # first read the taxonomy mapping
        open( TAXONMAP, "$markerdir/$marker.ncbimap") || croak("Unable to read file $markerdir/$marker.ncbimap\n");
	my %markerncbimap;
	while( my $line = <TAXONMAP> ){
                chomp($line);
		my ($markerbranch,$ncbiname) = split(/\t/, $line);
		$markerncbimap{$markerbranch} = [] unless defined( $markerncbimap{$markerbranch} );
		push( @{ $markerncbimap{$markerbranch} }, $ncbiname );
	}

        # then read & map the placement
        open(PLACEFILE, $placeFile) || croak("Unable to read file $placeFile\n");
	my $placeline = 0;
	my %curplaces;
	while( my $line = <PLACEFILE> ){
		$placeline=1 if($line =~ /"placements"/);
		next if($line =~ /^\>/);
		next if($line =~ /^\s*\#/);
		if($placeline==1 && $line =~ /\[(\d+),\s.\d+\.?\d+,\s(\d+\.?\d*),/){
			$curplaces{$1}=$2;
		}
		if($placeline==1 && $line =~ /\"n\"\:\s+\[\"(.+?)\"\]/){
			my $qname = $1;
			foreach my $edge(keys(%curplaces)){
				my $weightRatio = $curplaces{$edge};
				$weightRatio *= $coverage{$qname} if defined($coverage{$qname});
				my $mapcount = scalar(@{$markerncbimap{$edge}});
				foreach my $taxon( @{$markerncbimap{$edge}} ){
					my ($taxon_name, $taxon_level, $taxon_id) = getTaxonInfo($taxon);
					$placements{$qname} = () unless defined($placements{$qname});
					$placements{$qname}{$taxon_id} = 0 unless defined($placements{$qname}{$taxon_id});
					$placements{$qname}{$taxon_id} += $curplaces{$edge} / $mapcount;

					$ncbireads{$taxon} = 0 unless defined $ncbireads{$taxon};
					$ncbireads{$taxon} += $weightRatio / $mapcount;  # split the p.p. across the possible edge mappings
				}
			}
			%curplaces = ();
		}
	}
}
    # also write out the taxon assignments for sequences
    open(SEQUENCETAXA, ">".$self->{"fileDir"}."/sequence_taxa.txt");
    open(SEQUENCESUMMARY, ">".$self->{"fileDir"}."/sequence_taxa_summary.txt");
	foreach my $qname( keys( %placements ) ){
		# sum up all placements for this sequence, use to normalize
		my $placecount = 0;
		foreach my $taxon_id( keys( %{$placements{$qname}} ) ){
			$placecount += $placements{$qname}{$taxon_id};
		}
		# normalize to probability distribution
		foreach my $taxon_id (sort {$placements{$qname}{$b} <=> $placements{$qname}{$a} } keys %{$placements{$qname}}) {
			$placements{$qname}{$taxon_id} /= $placecount;			
			my ($taxon_name, $taxon_level, $tid) = getTaxonInfo($taxon_id);
			print SEQUENCETAXA "$qname\t$taxon_id\t$taxon_level\t$taxon_name\t".$placements{$qname}{$taxon_id}."\n";
		}
		my $readsummary = sum_taxon_levels( $placements{$qname} );
		foreach my $taxon_id (sort {$readsummary->{$b} <=> $readsummary->{$a} } keys %{$readsummary}) {
			my ($taxon_name, $taxon_level, $tid) = getTaxonInfo($taxon_id);
			print SEQUENCESUMMARY "$qname\t$taxon_id\t$taxon_level\t$taxon_name\t".$readsummary->{$taxon_id}."\n";
		}
		
	}
    close(SEQUENCESUMMARY);
    close(SEQUENCETAXA);
	
    # sort descending
    open(taxaOUT,">".$self->{"fileDir"}."/taxasummary.txt");
    foreach my $taxon (sort {$ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads) {
	my ($taxon_name, $taxon_level, $taxon_id) = getTaxonInfo($taxon);
	print taxaOUT join("\t",$taxon_id,$taxon_level,$taxon_name, $ncbireads{$taxon}),"\n";
    }
    close(taxaOUT);
    
    # sample from multinomial to get confidence limits
    # get total read count
    my $totalreads=0;
    foreach my $val(values(%ncbireads)){
        $totalreads+=$val;
    }
    print STDERR "Total reads are $totalreads\n";

    # write the taxa with 90% highest posterior density, assuming each read is an independent observation
    my $taxasum = 0;
    open(TAXAHPDOUT,">".$self->{"fileDir"}."/taxa_90pct_HPD.txt");
    foreach my $taxon (sort {$ncbireads{$b} <=> $ncbireads{$a} } keys %ncbireads) {
	$taxasum += $ncbireads{$taxon};
	my ($taxon_name, $taxon_level, $taxon_id) = getTaxonInfo($taxon);
	print TAXAHPDOUT join("\t",$taxon_id,$taxon_level,$taxon_name, $ncbireads{$taxon}),"\n";
	last if $taxasum >= $totalreads*0.9;
    }
    close(TAXAHPDOUT);

}

#
# non-functional until dependency on Math::Random can be eliminated
#
sub write_confidence_intervals{
    my $self = shift;
    my $ncbireadsref = shift;
    my $totalreads = shift;
    my %ncbireads = %$ncbireadsref;
    # normalize to a sampling distribution
    foreach my $key(keys(%ncbireads)){
        $ncbireads{$key}/=$totalreads + 1;
    }
    my $normsum=0;
    my @valarray = values(%ncbireads);
    foreach my $val(@valarray){
        $normsum+=$val;
    }

    my $sample_count = 100;
    my %samples;
    for( my $sI=0; $sI<$sample_count; $sI++ ){
#        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
	my @sample;
        my $kI=0;
        foreach my $key(keys(%ncbireads)){
            push(@{ $samples{$key} }, $sample[$kI++]);
        }
    }
    
    open(taxaCONF,">".$self->{"fileDir"}."/taxaconfidence.txt");
    foreach my $key(keys(%samples)){
        my @svals = @{ $samples{$key} };
        my @sorted = sort { $a <=> $b } @svals;
	my ($taxon_name, $taxon_level, $taxon_id) = getTaxonInfo($key);
        print taxaCONF join("\t", $taxon_id, $taxon_level, $taxon_name, $sorted[0], $sorted[int($sample_count*0.1)], $sorted[int($sample_count*0.25)], $sorted[int($sample_count*0.5)], $sorted[int($sample_count*0.75)], $sorted[int($sample_count*0.9)], $sorted[$sample_count-1]), "\n";
    }
}

sub sum_taxon_levels{
	my $placements = shift;
	my %summarized = ();
	foreach my $taxon_id( keys %$placements ){
		my $cur_tid = $taxon_id;		
		while( defined($cur_tid) && $cur_tid != 1 ){
			$summarized{$cur_tid} = 0 unless defined($summarized{$cur_tid});			
			$summarized{$cur_tid} += $placements->{$taxon_id};
			$cur_tid = $parent{$cur_tid}[0];
		}
	}
	return \%summarized;
}

sub getTaxonInfo {
	my $in = shift;
	if($in =~ /^\d+$/){
		#it's an ncbi taxon id.  look up its name and level.
		my $name = $idnamemap{$in};
		my $level = $parent{$in}->[1];
		return ($name, $level, $in);
	}elsif($in =~ /\w+/){
		# old style map, need to go from NCBI name back to ID
		my $name = $in;
		$name =~ s/_/ /g; # map uses spaces instead of underscores
		my ($id, $qname) = dongyingFindNameInTaxaDb($name);
		my $level = $parent{$id}->[1];
		return ($in, $level, $id);
	}else{
	}
	return ($in,"","");
}

sub treeName {
    my $inName = shift;
    $inName=~s/\s+/_/g;
    $inName=~s/'//g;
    $inName=~s/[\(\)]//g;
    $inName=~s/-/_/g;
    $inName=~s/\//_/g;
    return $inName;
}

=head2 homogenizeNameAlaDongying

=cut

sub homogenizeNameAlaDongying {
    my $inName = shift;
    return "" unless defined($inName);
    $inName=~s/^\s+//;
    $inName=~s/\s+$//;
    $inName=~s/\s+/ /g;
    $inName=~s/,//g;
    $inName=~s/[)(]//g;
    $inName=~s/-/ /g;
    $inName=uc $inName;
    return $inName;
}

=head2 dongyingFindNameInTaxaDb
    
=cut
    
sub dongyingFindNameInTaxaDb {
    my $name = shift;
    return "" unless defined($name);
    $name=~s/^\s+//;
    my @t=split(/\s+/, $name);
    my $input_name=join(" ",@t);
    my $q_name=$input_name;
    my $id="ERROR";
    while(@t>=1){
	$q_name=join(" ",@t);
	$q_name=uc $q_name;
	if(defined($nameidmap{$q_name})){
	    $id=$nameidmap{$q_name};
	    last;
	}
	pop @t;
    }
    return ($id,$q_name);
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

    perldoc Amphora2::Summarize


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

1; # End of Amphora2::Summarize.pm

