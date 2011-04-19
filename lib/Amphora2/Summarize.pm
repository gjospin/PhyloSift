package Amphora2::Summarize;

use warnings;
use strict;
use autodie;

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

sub summarize {
    @ARGV = @_;
    my $markerdir = "markers";
    open( my $NAMETABLE, "$markerdir/name.table" );
    my %namemap;
    while( my $line = <$NAMETABLE> ){
	chomp $line;
	my @vals = split( /\t/, $line );
	$namemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]);
    }
    
    open( my $TAXIDS, "ncbi/names.dmp" );
    while( my $line = <$TAXIDS> ){
	chomp $line;
	if(($line =~ /scientific name/) || ($line =~ /synonym/) || ($line =~ /misspelling/)){
	    my @vals = split( /\s+\|\s+/, $line );
	    $nameidmap{homogenizeNameAlaDongying($vals[1])}=$vals[0];
	    $idnamemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]) if($line =~ /scientific name/);
	}
    }
    
    open( my $TAXSTRUCTURE, "ncbi/nodes.dmp" );
    my %parent;
    while( my $line = <$TAXSTRUCTURE> ){
	chomp $line;
	my @vals = split( /\s+\|\s+/, $line );
	$parent{$vals[0]} = [$vals[1],$vals[2]];
    }
    
    my %hitcounter;
    my $readcount = 0;
    open(my $neighborIN, $ARGV[0]);
    
    while( my $line = <$neighborIN> ){
	chomp $line;
	$line =~ s/\s+$//g;
	$line =~ s/^\s+//g;
	#	next unless length($line) > 12;
	next unless( defined($namemap{$line}) );
	my ($tid,$name) = dongyingFindNameInTaxaDb($namemap{$line});
	if($tid eq "ERROR"){
	    print STDERR "Error! Could not find $line in name map\n" if length($line) > 12;
	    next;
	}
	
		#got the taxon id, now walk to root tallying everything we hit
	next unless(defined($tid));
	while( $tid != 1 ){
	    if(defined($hitcounter{$tid})){
		$hitcounter{$tid}++;
	    }else{
		$hitcounter{$tid}=1;
	    }
	    $tid = $parent{$tid}->[0];
	}
	$readcount++;
    }
    close($neighborIN);
    my %hitvals;
    foreach my $tid(keys(%hitcounter)){
	my $frac = sprintf("%.4f",$hitcounter{$tid}/$readcount);
	#	$hitvals{$idnamemap{$tid}} = $frac;
	$hitvals{$idnamemap{$tid}}=[$hitcounter{$tid},$parent{$tid}->[1],$frac];
	
    }
    my @sorted = reverse sort { $hitvals{$a}->[0] <=> $hitvals{$b}->[0] } keys %hitvals; 
    open(my $taxaOUT,">$ARGV[1]");
    foreach my $names (@sorted){
	print taxaOUT join("\t",$hitvals{$names}->[1],$names,$hitvals{$names}->[0],$hitvals{$names}->[2]),"\n";
    }
    close($taxaOUT);
}

=head2 dongyingFindNameInTaxaDb

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

