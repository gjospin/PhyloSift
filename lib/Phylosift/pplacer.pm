package Phylosift::pplacer;
use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(debug);
use Phylosift::Summarize;
use Bio::Phylo::IO qw(parse unparse);

=head1 NAME

Phylosift::pplacer - place aligned reads onto a phylogenetic tree with pplacer

=head1 VERSION

Version 0.01

=cut
our $VERSION = '0.01';

=head1 SYNOPSIS

Runs Pplacer for all the markers in the list passed as input.
the script will look for the files of interest in
$workingDir/markers/
AND
$workingDir/PS_temp/alignDir/

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 pplacer

=cut

sub pplacer {
	my %args    = @_;
	my $self    = $args{self} // miss("self");
	my $markRef = $args{marker_reference} // miss("marker_reference");
	directoryPrepAndClean( self => $self );

	# if we have a coverage map then weight the placements
	my $covref;
	if ( defined( $self->{"coverage"} ) ) {
		$covref = Phylosift::Summarize::read_coverage( file => $self->{"coverage"} );
	}
	if ( $self->{"updated"} ) {
		my $markerPackage = Phylosift::Utilities::get_marker_package( self => $self, marker => "concat" );
		my $pp =
		    "$Phylosift::Utilities::pplacer --verbosity 0 -c $markerPackage -j "
		  . $self->{"threads"}
		  . " --groups 5 "
		  . $self->{"alignDir"}
		  . "/concat.trim.fasta";
		system($pp);
		`mv $self->{"workingDir"}/concat.trim.jplace $self->{"treeDir"}` if ( -e $self->{"workingDir"} . "/concat.trim.jplace" );
		return unless -e $self->{"treeDir"} . "/concat.trim.jplace";
		name_taxa_in_jplace( self => $self, input => $self->{"treeDir"} . "/concat.trim.jplace", output => $self->{"treeDir"} . "/concat.trim.jplace" );
		return unless defined($covref);
		weight_placements( self => $self, coverage => $covref, place_file => $self->{"treeDir"} . "/concat.trim.jplace" );
		return;
	}
	foreach my $marker ( @{$markRef} ) {
		my $readAlignmentFile = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker );
		my $readAlignmentDNAFile = $self->{"alignDir"} . "/" . Phylosift::Utilities::get_aligner_output_fasta_DNA( marker => $marker );
		next unless -e $readAlignmentFile || -e $readAlignmentDNAFile;
		my $markerPackage = Phylosift::Utilities::get_marker_package( self => $self, marker => $marker );
		debug "Running Placer on $marker ....\t";
		my $placeFile = Phylosift::Utilities::get_read_placement_file( marker => $marker );
		my $placeFileDNA = Phylosift::Utilities::get_read_placement_file_DNA( marker => $marker );
		if ( $self->{"updated"} == 0 ) {
			my $pp = "";
			if ( Phylosift::Utilities::marker_oldstyle( markers => $marker ) ) {

				# run pplacer the old way, using phyml trees which aren't supported by reference packages
				my $trimfinalFastaFile = Phylosift::Utilities::get_trimfinal_fasta_marker_file( self => $self, marker => $marker );
				my $trimfinalFile = Phylosift::Utilities::get_trimfinal_marker_file( self => $self, marker => $marker );
				my $treeFile = Phylosift::Utilities::get_tree_marker_file( self => $self, marker => $marker );
				my $treeStatsFile = Phylosift::Utilities::get_tree_stats_marker_file( self => $self, marker => $marker );

				# Pplacer requires the alignment files to have a .fasta extension
				if ( !-e "$trimfinalFastaFile" ) {
					`cp $trimfinalFile $trimfinalFastaFile`;
				}
				$pp =
				    "$Phylosift::Utilities::pplacer --verbosity 0 -j "
				  . $self->{"threads"}
				  . " -r $trimfinalFastaFile -t $treeFile -s $treeStatsFile $readAlignmentFile";
			} else {
				$pp = "$Phylosift::Utilities::pplacer --verbosity 0 -c $markerPackage -j " . $self->{"threads"} . " $readAlignmentFile";
			}
			debug "Running $pp\n";
			system("$pp");
		} else {

			#run pplacer on amino acid data
			my $pp = "$Phylosift::Utilities::pplacer --verbosity 0 -j " . $self->{"threads"} . " -c $markerPackage $readAlignmentFile";
			print "Running $pp\n";
			system($pp);

			#run pplacer on nucleotide data
			if ( -e $readAlignmentDNAFile ) {
				my $codonmarkers = $markerPackage;
				$codonmarkers =~ s/.updated/.codon.updated/g;
				my $pp = "$Phylosift::Utilities::pplacer --verbosity 0 -j " . $self->{"threads"} . " -c $codonmarkers $readAlignmentDNAFile";
				print "Running $pp\n";
				system($pp);
			}
		}

		# pplacer writes its output to the directory it is called from. Need to move the output to the trees directory
		`mv $self->{"workingDir"}/$placeFile $self->{"treeDir"}`    if ( -e $self->{"workingDir"} . "/$placeFile" );
		`mv $self->{"workingDir"}/$placeFileDNA $self->{"treeDir"}` if ( -e $self->{"workingDir"} . "/$placeFileDNA" );
		name_taxa_in_jplace( self => $self, input => $self->{"treeDir"} . "/$placeFile", output => $self->{"treeDir"} . "/$placeFile" )
		  if ( -e $self->{"treeDir"} . "/$placeFile" );
		name_taxa_in_jplace( self => $self, input => $self->{"treeDir"} . "/$placeFileDNA", output => $self->{"treeDir"} . "/$placeFileDNA" )
		  if ( -e $self->{"treeDir"} . "/$placeFileDNA" );
		next unless ( defined($covref) );
		weight_placements( self => $self, coverage => $covref, place_file => $self->{"treeDir"} . "/$placeFile" );
	}
}

=head2 weight_placements

=cut

sub weight_placements {
	my %args       = @_;
	my $coverage   = $args{coverage} // miss("coverage");
	my $place_file = $args{place_file} // miss("place_file");

	# weight the placements
	my $INPLACE = ps_open(  $place_file );
	my $OUTPLACE = ps_open( ">$place_file.wt" );
	my $placeline = 0;
	while ( my $line = <$INPLACE> ) {
		$placeline = 1 if ( $line =~ /"placements"/ );

		# have we reached the end of a placement entry?
		if ( $placeline == 1 && $line =~ /\"n\"\:\s+\[\"(.+?)\"\]/ ) {
			my $qname = $1;
			print STDERR "Found placeline for $qname\n";
			my @qnames = split( /,/, $qname );

			# if we have read coverage information, add it to an updated placement file
			my $mass = 0;
			foreach my $n ($qnames) {
				$mass += $coverage->{$n} if defined( $coverage->{$n} );
				print STDERR "Unable to find coverage for $qname\n" unless defined( $coverage->{$n} );
			}
			print $OUTPLACE ", \"m\": \"$mass\",\n";
		}
		if ( $placeline == 1 && $line =~ /\"nm\"\:\s+\[\[\"(.+?)\",/ ) {
			my $qname = $1;
			$qname =~ s/\\\/\d-\d+//g;
			print STDERR "Found placeline for $qname\n";

			# if we have read coverage information, add it to an updated placement file
			my $mass = 0;
			$mass += $coverage->{$qname} if defined( $coverage->{$qname} );
			$line =~ s/\"nm\"\:\s+\[\[\"(.+?)\", \d+/"nm": [[\"$qname\", $mass/g;
			print STDERR "Unable to find coverage for $qname\n" unless defined( $coverage->{$qname} );
		}
		print $OUTPLACE $line;
	}
	close($OUTPLACE);
	`mv $place_file.wt $place_file`;
}

=head2 directoryPrepAndClean

=cut

sub directoryPrepAndClean {
	my %args = @_;
	my $self = $args{self} // miss("self");
	`mkdir $self->{"tempDir"}` unless ( -e $self->{"tempDir"} );

	#create a directory for the Reads file being processed.
	`mkdir $self->{"fileDir"}` unless ( -e $self->{"fileDir"} );
	`mkdir $self->{"treeDir"}` unless ( -e $self->{"treeDir"} );
}

=head1 SUBROUTINES/METHODS

=head2 nameTaxa

=cut

sub name_taxa_in_jplace {
	my %args   = @_;
	my $self   = $args{self} // miss("self");
	my $input  = $args{input} // miss("input");
	my $output = $args{output} // miss("output");

	# read in the taxon name map
	my %namemap;
	my $taxon_map_file = Phylosift::Utilities::get_marker_taxon_map( self => $self );
	my $NAMETABLE = ps_open( $taxon_map_file );
	while ( my $line = <$NAMETABLE> ) {
		chomp $line;
		my @pair = split( /\t/, $line );
		$namemap{ $pair[0] } = $pair[1];
	}
	Phylosift::Summarize::read_ncbi_taxon_name_map();

	# parse the tree file to get leaf node names
	# replace leaf node names with taxon labels
	my $TREEFILE = ps_open( $input );
	my @treedata = <$TREEFILE>;
	close $TREEFILE;
	my $tree_string = $treedata[1];
	$tree_string =~ s/^\s+\"//g;
	$tree_string =~ s/\{\d+?\}//g;
	$tree_string =~ s/\"\,$//g;
	my $tree = Bio::Phylo::IO->parse(
									  '-string' => $tree_string,
									  '-format' => 'newick',
	)->first;

	foreach my $node ( @{ $tree->get_entities } ) {

		# skip this one if it is not a leaf
		#		next if ( scalar($node->get_children())>0 );
		my $name = $node->get_name;
		next unless defined $namemap{$name};
		my @data = Phylosift::Summarize::get_taxon_info( taxon => $namemap{$name} );
		my $ncbi_name = Phylosift::Summarize::tree_name( name => $data[0] );
		$node->set_name($ncbi_name);
	}
	my $new_string = "  \"" . unparse( '-phylo' => $tree, '-format' => 'newick' ) . "\",\n";
	$treedata[1] = $new_string;
	$TREEFILE = ps_open( ">$output" );
	print $TREEFILE @treedata;
	close $TREEFILE;
}

=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::pplacer


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Phylosift-Phylosift>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Phylosift-Phylosift>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Phylosift-Phylosift>

=item * Search CPAN

L<http://search.cpan.org/dist/Phylosift-Phylosift/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation.

See http://dev.perl.org/licenses/ for more information.


=cut
1;    # End of Phylosift::pplacer.pm
