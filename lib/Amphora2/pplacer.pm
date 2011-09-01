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
    #debug "PPLACER MARKS\t @{$markRef}\n";


    foreach my $marker(@{$markRef}){

	my $trimfinalFastaFile = Amphora2::Utilities::getTrimfinalFastaMarkerFile($self,$marker);
	my $trimfinalFile = Amphora2::Utilities::getTrimfinalMarkerFile($self,$marker);
	my $treeFile = Amphora2::Utilities::getTreeMarkerFile($self,$marker);
	my $treeStatsFile = Amphora2::Utilities::getTreeStatsMarkerFile($self,$marker);
	
	# Pplacer requires the alignment files to have a .fasta extension
	if(!-e "$Amphora2::Utilities::marker_dir/$trimfinalFastaFile" ){
	    `cp $Amphora2::Utilities::marker_dir/$trimfinalFile $Amphora2::Utilities::marker_dir/$trimfinalFastaFile`;
	}
#	if(!-e $self->{"alignDir"}."/$marker.aln_hmmer3.trim.fasta"  && -e $self->{"alignDir"}."/$marker.aln_hmmer3.trim"){
#	    `cp $self->{"alignDir"}/$marker.aln_hmmer3.trim $self->{"alignDir"}/$marker.aln_hmmer3.trim.fasta`;
#	}else{
#	    next;
#	}
	#adding a printed statement to check on the progress
	print STDERR "Running Placer on $marker ....\t";
	#running Pplacer
	if(!-e $self->{"treeDir"}."/$marker.aln_hmmer3.trim.place"){
	    `$Amphora2::Utilities::pplacer -p -r $Amphora2::Utilities::marker_dir/$trimfinalFastaFile -t $Amphora2::Utilities::marker_dir/$treeFile -s $Amphora2::Utilities::marker_dir/$treeStatsFile $self->{"alignDir"}/$marker.aln_hmmer3.trim`;
	}
	#adding a printed statement to check on the progress (not really working if using parrallel jobs)
	print STDERR "Done !\n";
	#Pplacer write its output to the directory it is called from. Need to move the output to the trees directory
	if(-e $self->{"workingDir"}."/$marker.aln_hmmer3.trim.place"){
	    `mv $self->{"workingDir"}/$marker.aln_hmmer3.trim.place $self->{"treeDir"}`;
	}
	
	#Transform the .place file into a tree file
	if(-e $self->{"treeDir"}."/$marker.aln_hmmer3.trim.place" && !-e $self->{"treeDir"}."/$marker.aln_hmmer3.trim.tog.tree"){
	    `placeviz --loc -p $self->{"treeDir"}/$marker.aln_hmmer3.trim.place`
	}
	
	
	#placeviz writes its output to the directory it was called from, need to move the output to the trees directory
	my $workingDir = $self->{"workingDir"};
	my @placevizFiles = <$workingDir/$marker.aln_hmmer3.*>;
	if(scalar(@placevizFiles)>0){
	    `mv $self->{"workingDir"}/$marker.* $self->{"treeDir"}`;
	}
	#added the .PP. check to accomodate for what pplacer names its files (tax branch or master branch)
	# transform the taxon names in the tree file
	if(-e $self->{"treeDir"}."/$marker.aln_hmmer3.trim.tog.tre"){
	    nameTaxa($self->{"treeDir"}."/$marker.aln_hmmer3.trim.tog.tre");
	}elsif(-e $self->{"treeDir"}."/$marker.aln_hmmer.trim.PP.tog.tre"){
	    nameTaxa($self->{"treeDir"}."/$marker.aln_hmmer3.trim.PP.tog.tre");
	}
	
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
