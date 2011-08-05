package Amphora2::Summarize;

use warnings;
use strict;
use Amphora2::Amphora2;
use Math::Random;
use Carp;

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
    open( my $TAXIDS, "$ncbidir/names.dmp" );
    while( my $line = <$TAXIDS> ){
	chomp $line;
	if(($line =~ /scientific name/) || ($line =~ /synonym/) || ($line =~ /misspelling/)){
	    my @vals = split( /\s+\|\s+/, $line );
	    $nameidmap{homogenizeNameAlaDongying($vals[1])}=$vals[0];
	    $idnamemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]) if($line =~ /scientific name/);
	}
    }
}

# now read the NCBI taxonomy structure
# puts the results in a hash called "parent"
my %parent;
sub readNcbiTaxonomyStructure {
    my $ncbidir = $Amphora2::Utilities::ncbi_dir;
    open( my $TAXSTRUCTURE, "$ncbidir/nodes.dmp" );
    while( my $line = <$TAXSTRUCTURE> ){
	chomp $line;
	my @vals = split( /\s+\|\s+/, $line );
	$parent{$vals[0]} = [$vals[1],$vals[2]];
    }
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
	
    # read all of the .place files for markers
    # map them onto the ncbi taxonomy
    foreach my $marker(@{$markRef}){
	# don't bother with this one if there's no read placements
	next unless( -e $self->{"treeDir"}."/$marker.aln_hmmer3.place" );

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
        open(PLACEFILE, $self->{"treeDir"}."/$marker.aln_hmmer3.place") || croak("Unable to read file ".$self->{"treeDir"}."/$marker.aln_hmmer3.place\n");
	my $placeline = 0;
	while( my $line = <PLACEFILE> ){
            $placeline=1 if($line =~ /^\>/);
            next if($line =~ /^\>/);
            next if($line =~ /^\s*\#/);
            next unless($line =~ /^\d+\t\d+/);
            if($placeline==1){
                my @pline = split(/\t/, $line);
                my $mapcount = scalar(@{$markerncbimap{$pline[0]}});
                foreach my $taxon( @{$markerncbimap{$pline[0]}} ){
                    $ncbireads{$taxon} = 0 unless defined $ncbireads{$taxon};
                    $ncbireads{$taxon} += $pline[1] / $mapcount;	# split the p.p. across the possible edge mappings
                }
            }
	}
    }
    open(taxaOUT,">".$self->{"fileDir"}."/taxasummary.txt");
    foreach my $taxon(keys(%ncbireads)){
	print taxaOUT join("\t",$taxon,$ncbireads{$taxon}),"\n";
    }
    close(taxaOUT);
    
    # sample from multinomial to get confidence limits
    # get total read count
    my $totalreads=0;
    foreach my $val(values(%ncbireads)){
        $totalreads+=$val;
    }
    print "Have $totalreads reads\n";
    # normalize to a sampling distribution
    foreach my $key(keys(%ncbireads)){
        $ncbireads{$key}/=$totalreads + 1;
    }
    my $normsum=0;
    my @valarray = values(%ncbireads);
    foreach my $val(@valarray){
        $normsum+=$val;
    }
    print "Normalize sum $normsum reads\n";
#    $valarray[0] += 1 - $normsum; # ugh, deal with fp error
    my $sample_count = 100;
    my %samples;
    for( my $sI=0; $sI<$sample_count; $sI++ ){
        my @sample = Math::Random::random_multinomial( $totalreads, @valarray );
        my $kI=0;
        foreach my $key(keys(%ncbireads)){
            push(@{ $samples{$key} }, $sample[$kI++]);
        }
    }
    
    open(taxaCONF,">".$self->{"fileDir"}."/taxaconfidence.txt");
    foreach my $key(keys(%samples)){
        my @svals = @{ $samples{$key} };
        my @sorted = sort { $a <=> $b } @svals;
        print taxaCONF join("\t",$key, $sorted[0], $sorted[int($sample_count*0.1)], $sorted[int($sample_count*0.25)], $sorted[int($sample_count*0.5)], $sorted[int($sample_count*0.75)], $sorted[int($sample_count*0.9)], $sorted[$sample_count-1]), "\n";
    }
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
    $inName=~s/^\s+//;
    $inName=~s/\s+$//;
    $inName=~s/\s+/ /g;
    $inName=~s/,//g;
    $inName=uc $inName;
    return $inName;
}

=head2 dongyingFindNameInTaxaDb
    
=cut
    
sub dongyingFindNameInTaxaDb {
    my $name = shift;
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

