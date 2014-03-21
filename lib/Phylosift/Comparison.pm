package Phylosift::Comparison;
use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Phylosift::Phylosift;
use Phylosift::Utilities;
use Bio::Phylo::IO qw(parse unparse);

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::Comparison - compares phylogenetic structure across many samples

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS

Given a set of samples, will use guppy to calculate Kantorovich-Rubenstein distance among samples
and also do an edge PCA on them. Transforms the XML output into something with meaningful taxon names.

=head1 SUBROUTINES/METHODS

=head2 compare

=cut

sub compare {
	my %args             = @_;
	my $self             = $args{self} || miss("self");
	my $parent_directory = $args{parent_dir} || miss("parent_dir");

	# what do we want to accomplish with this?
	# simplest approach:
	# 1. take many jplace files on concat alignments
	# 2a. compute squash clustering
	# 2b. calculate edge PCA and create tree XML
	# 3. swap real names back into tree XML
	# 4. make some pretty pictures of tree, PCA plot, and squash cluster
	#
	# more complicated, possibly better:
	# 1. merge all jplaces to a single file
	# 2. use rppr voronoi to select a small set of taxa that are relevant for this set of samples
	# 3. create a new concat package with the rppr voronoi taxa
	# 4. re-place reads from each sample
	# 5. now do the above steps listed in "simple approach"
	my $JPLACES = ps_open("ls -1 $parent_directory/*/treeDir/concat.trim.jplace |");
	my @files;
	while ( my $file = <$JPLACES> ) {
		chomp $file;
		push( @files, $file );
	}
	my $squash_cl = "$Phylosift::Utilities::guppy squash --prefix squash ".join( " ", @files );
	print STDERR "Squashing with: $squash_cl ";
	system($squash_cl);

	#	output goes to cluster.tre
	my $pca_cl = "$Phylosift::Utilities::guppy pca --prefix pca ".join( " ", @files );
	print STDERR "pca with: $pca_cl ";
	system($pca_cl);
	my $merge_cl = "$Phylosift::Utilities::guppy merge -o $parent_directory/merged.jplace ".join( " ", @files );
	system($merge_cl);
	my @processfiles = @files;
	unshift( @files, "$parent_directory/merged.jplace" );
	foreach my $file (@files) {
		my $rarefact_cl = "$Phylosift::Utilities::guppy rarefact -o \"$file.rarefaction\" \"$file\"";
		system($rarefact_cl);
		my $compress_cl = "$Phylosift::Utilities::guppy compress -o \"$file.compression\" --cutoff 0.4 \"$file\"";
		system($compress_cl);
	}
	rename_nodes( in_file => "pca.xml", out_file => "pca.named.xml" );
}

sub rename_nodes {
	my %args    = @_;
	my $infile  = $args{in_file} || miss("in_file");
	my $outfile = $args{out_file} || miss("out_file");
	my $data = parse(
					  '-file'       => $infile,
					  '-format'     => 'phyloxml',
					  '-as_project' => 1
	);
	my ($taxa) = @{ $data->get_taxa };
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

    perldoc Phylosift::Comparison


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

1;    # End of Phylosift::Comparison.pm
