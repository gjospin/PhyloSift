package Amphora2::pplacer;

use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Amphora2::Amphora2;
=head1 NAME

Amphora2::pplacer - place aligned reads onto a phylogenetic tree with pplacer

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Runs Pplacer for all the markers in the list passed as input.
the script will look for the files of interest in
$workingDir/markers/
AND
$workingDir/Amph_temp/alignments/

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 pplacer

=cut

sub pplacer {
    my $self = shift;
    my $markRef = shift;
    directoryPrepAndClean($self);


    foreach my $marker(@{$markRef}){

	my $readAlignmentFile = $self->{"alignDir"}."/".Amphora2::Utilities::getAlignerOutputFastaAA($marker);
	my $readAlignmentDNAFile = $self->{"alignDir"}."/".Amphora2::Utilities::getAlignerOutputFastaDNA($marker);
	my $markerPackage = Amphora2::Utilities::getMarkerPackage($self,$marker);
	
	print STDERR "Running Placer on $marker ....\t";
	my $placeFile = Amphora2::Utilities::getReadPlacementFile($marker);
	my $placeFileDNA = Amphora2::Utilities::getReadPlacementFileDNA($marker);
	if($self->{"updated"}==0){
		# run pplacer the old way, using phyml trees which aren't supported by reference packages
		my $trimfinalFastaFile = "$Amphora2::Utilities::marker_dir/".Amphora2::Utilities::getTrimfinalFastaMarkerFile($self,$marker);
		my $trimfinalFile = "$Amphora2::Utilities::marker_dir/".Amphora2::Utilities::getTrimfinalMarkerFile($self,$marker);
		my $treeFile = "$Amphora2::Utilities::marker_dir/".Amphora2::Utilities::getTreeMarkerFile($self,$marker);
		my $treeStatsFile = "$Amphora2::Utilities::marker_dir/".Amphora2::Utilities::getTreeStatsMarkerFile($self,$marker);
		
		# Pplacer requires the alignment files to have a .fasta extension
		if(!-e "$trimfinalFastaFile" ){
		    `cp $trimfinalFile $trimfinalFastaFile`;
		}
	
	    my $pp = "$Amphora2::Utilities::pplacer -p -r $trimfinalFastaFile -t $treeFile -s $treeStatsFile $readAlignmentFile";
	    `$pp`;
	}else{
		#run pplacer on amino acid data
		my $pp = "$Amphora2::Utilities::pplacer -c $markerPackage --no-pre-mask $readAlignmentFile";
		print "Running $pp\n";
	    system($pp);
		#run pplacer on nucleotide data
		if(-e $readAlignmentDNAFile){
			my $codonmarkers = $markerPackage;
			$codonmarkers =~ s/.updated/.codon.updated/g;
			my $pp = "$Amphora2::Utilities::pplacer -c $codonmarkers --no-pre-mask $readAlignmentDNAFile";
			print "Running $pp\n";
		    system($pp);
		}
	}
	print STDERR "Done !\n";
	# pplacer writes its output to the directory it is called from. Need to move the output to the trees directory
    `mv $self->{"workingDir"}/$placeFile $self->{"treeDir"}` if(-e $self->{"workingDir"}."/$placeFile");
    `mv $self->{"workingDir"}/$placeFileDNA $self->{"treeDir"}` if(-e $self->{"workingDir"}."/$placeFileDNA");
    }
}


=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean{

    my $self = shift;
    my @markers = @_;
    `mkdir $self->{"tempDir"}` unless (-e $self->{"tempDir"});
    #create a directory for the Reads file being processed.
    `mkdir $self->{"fileDir"}` unless (-e $self->{"fileDir"});
    `mkdir $self->{"treeDir"}` unless (-e $self->{"treeDir"});
}

=head1 SUBROUTINES/METHODS

=head2 nameTaxa

=cut


sub nameTaxa {
    my $filename = shift;

    # read in the taxon name map
    my %namemap;
    open(NAMETABLE, "$Amphora2::Utilities::marker_dir/name.table")or die "Couldn't open $Amphora2::Utilities::marker_dir/name.table\n";
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

=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-amphora2-amphora2 at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Amphora2-Amphora2>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Amphora2::pplacer


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

1; # End of Amphora2::pplacer.pm
