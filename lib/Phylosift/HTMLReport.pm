package Phylosift::HTMLReport;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Bio::AlignIO;
use Phylosift::Phylosift;
use Phylosift::Utilities qw(:all);
use Phylosift::Summarize;
use Phylosift::Settings;
use Bio::Phylo;
use Bio::Phylo::Forest::Tree;
use Bio::Phylo::IO qw(parse unparse);
use Cwd 'abs_path';
use Bio::AlignIO;
use XML::Writer;
use IO::File;
use File::Basename;
use Carp;

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::HTMLReport - support for generating HTML reports summarizing results

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS


=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 begin_report

=cut

sub begin_report {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $file   = $args{file} || "krona.html";
	my $OUTPUT = IO::File->new(">$file");
	print $OUTPUT <<EOF;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
 <head>
  <meta charset="utf-8"/>
  <base href="http://krona.sourceforge.net/" target="_blank"/>
  <link rel="shortcut icon" href="img/favicon.ico"/>
  <script id="notfound">window.onload=function(){document.body.innerHTML="Could not get resources from \"http://krona.sourceforge.net\"."}</script>
  <script src="src/krona-2.0.js"></script>
 </head>
 <body>
  <img id="hiddenImage" src="img/hidden.png" style="display:none"/>
  <noscript>Javascript must be enabled to view this page.</noscript>
  <div style=\"position:absolute;bottom:0;\">\n

EOF
	return $OUTPUT;
}

sub add_jnlp {
	my %args               = @_;
	my $self               = $args{self} || miss("self");
	my $marker             = $args{marker} || miss("marker");
	my $OUTPUT             = $args{OUTPUT} || miss("OUTPUT");
	my $xml                = $args{xml} || miss("xml");
	my $skip_jnlp_writeout = $args{html_only};
	$marker .= '.';
	$marker = '' if $marker eq "concat";
	my $jnlp     = $Phylosift::Settings::file_dir."/".$self->{"fileName"}.".$marker"."jnlp";
	my $xml_name = basename($xml);

	print $OUTPUT <<EOF;
<script>
var loc = window.location.pathname;
var dir = loc.substring(0, loc.lastIndexOf('/') + 1);
document.write('View phylogenetic placements for <a href="http://edhar.genomecenter.ucdavis.edu/cgi-bin/forester.jnlp?url='+dir+'$xml_name">marker $marker</a><br>');
</script>
EOF
	write_jnlp( self => $self, marker => $marker, jnlp => $jnlp, xml => $xml ) unless $skip_jnlp_writeout;
}

sub write_jnlp {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $marker = $args{marker} || miss("marker");
	my $xml    = $args{xml} || miss("xml");
	my $jnlp   = $args{jnlp};
	my $JOUT   = ps_open(">$jnlp");
	$xml = basename($xml);
	print $JOUT <<EOF;
<?xml version="1.0" encoding="UTF-8"?>
<jnlp spec="1.5+" codebase="http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/forester/">
    <information>
        <title>Archaeopteryx tree viewer</title>
        <vendor>Christian Zmasek, repackaged by Aaron Darling</vendor>
        <offline-allowed/>
    </information>
    <security>
      <all-permissions/>
    </security>
    <resources>
        <j2se version="1.5+" href="http://java.sun.com/products/autodl/j2se" java-vm-args="-Xmx300m"/>
        <jar href="http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/forester/forester.jar" main="true" />
    </resources>
    <application-desc main-class="org.forester.archaeopteryx.Archaeopteryx">
		<argument>$xml</argument>
    </application-desc>
    <update check="background"/>
</jnlp>	
EOF
	close($JOUT);
}

my $xml;
my %ncbi_summary;

sub add_krona {
	my %args            = @_;
	my $self            = $args{self} || miss("self");
	my $OUTPUT          = $args{OUTPUT} || miss("output");
	my $summary         = $args{summary} || miss("summary");
	my $KRONA_THRESHOLD = $Phylosift::Settings::krona_threshold;
	%ncbi_summary = ();
	%ncbi_summary = %$summary;
	print $OUTPUT "</div><div style=\"display:none\">\n<br>";
	$xml = new XML::Writer( OUTPUT => $OUTPUT );
	$xml->startTag( "krona", "collapse" => "false", "key" => "true" );
	$xml->startTag( "attributes", "magnitude" => "abundance" );
	$xml->startTag( "attribute",  "display"   => "reads" );
	$xml->characters("abundance");
	$xml->endTag("attribute");
	$xml->endTag("attributes");
	$xml->startTag("datasets");
	$xml->startTag("dataset");
	$xml->characters( $self->{"readsFile"} );
	$xml->endTag("dataset");
	$xml->endTag("datasets");
	debug "parse ncbi\n";

	# FIXME: work with other taxonomy trees
	my $taxonomy = Bio::Phylo::IO->parse(
										  '-file'   => "$Phylosift::Settings::marker_dir/ncbi_tree.updated.tre",
										  '-format' => 'newick',
	)->first;

	debug "visitor\n";
	my $root = $taxonomy->get_root;
	debug "Root node id ".$root->get_name."\n";
	debug "Root node read count ".$ncbi_summary{ $root->get_name }."\n";

	# write out abundance for nodes that have > $KRONA_THRESHOLD probability mass
	$root->visit_depth_first(
		-pre => sub {
			my $node = shift;
			my $name = $node->get_name;
			return
			  unless ( defined( $ncbi_summary{$name} )
					   && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD );
			my ( $taxon_name, $taxon_level, $taxon_id ) = Phylosift::Summarize::get_taxon_info( taxon => $name );
			$xml->startTag(
							"node",
							"name" => $taxon_name,
							"href" => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$name"
			);
			$xml->startTag("abundance");
			$xml->startTag("val");
			$xml->characters( $ncbi_summary{$name} );
			$xml->endTag("val");
			$xml->endTag("abundance");
		},
		-pre_sister => sub {
			my $node = shift;
			my $name = $node->get_name;
			return
			  unless ( defined( $ncbi_summary{$name} )
					   && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD );
			$xml->endTag("node");
		},
		-no_sister => sub {
			my $node = shift;
			my $name = $node->get_name;
			return
			  unless ( defined( $ncbi_summary{$name} )
					   && $ncbi_summary{$name} / $ncbi_summary{1} > $KRONA_THRESHOLD );
			$xml->endTag("node");
		}
	);
	debug "done visiting!\n";

	$xml->endTag("krona");
	$xml->end();
	print $OUTPUT "\n</div><br>";
}

sub add_run_info() {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $OUTPUT = $args{OUTPUT};
	print $OUTPUT "<div style=\"position:absolute;bottom:0;\">\n";
	Phylosift::Summarize::print_run_info( self => $self, OUTPUT => $OUTPUT, newline => "<br/>\n" );
	print $OUTPUT "</div></body></html>\n";
}

sub finalize {
	my %args   = @_;
	my $self   = $args{self} || miss("self");
	my $OUTPUT = $args{OUTPUT} || miss("OUTPUT");
	$OUTPUT->close();
}

=head1 AUTHOR

Aaron Darling, C<< <aarondarling at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::HTMLReport


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
