package Amphora2::Utilities;

use 5.006;
use strict;
use warnings;
use Carp;

=head1 NAME

Amphora2::Utilities - Implements miscellaneous accessory functions for Amphora2

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Amphora2::Amphora2;

    my $foo = Amphora2::Amphora2->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 programChecks

checks the program requirements for Amphora-2
writes to STDERR if a program is missing or the wrong version is installed
returns 1 or 0 depending on success of failure.

=cut

sub programChecks {
	eval 'require Bio::Seq;';
	if ($@) {
	    carp "Bioperl was NOT found\n";
	    return 1;
	}


	my $pplacercheck = `which pplacer`;
	if($pplacercheck eq ""){
	    #program not found return;
	    carp("pplacer v1.1.alpha00 not found");
	    return 1;
	}elsif(`pplacer -help` !~ m/v1.1.alpha00/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp("Warning : a different version of pplacer was found. Amphora-2 was tested with pplacer v1.1.alpha00\n");
	}else{
	    #program found and the correct version is installed
	}

	my $hmmercheck = `which hmmalign`;
	if($hmmercheck eq ""){
	    #program not found return;
	    carp("HMMER3 not found");
	    return 1;
	}elsif(`hmmalign -h` !~ m/HMMER 3.0rc1/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp "Warning : a different version of HMMER was found. Amphora-2 was tested with HMMER 3.0rc1\n";
	}else{
	    #program found and the correct version is installed
	}


	my $blastcheck = `which blastp`;
	if($blastcheck eq ""){
	    #program not found return;
	    carp("Blast 2.2.24+ not found");
	    return 1;
	}elsif(`blastp -help` !~ m/BLAST 2.2.24+/){
	    # pplacer was found but the version doens't match the one tested with Amphora
	    carp "Warning : a different version of Blast was found. Amphora-2 was tested with BLAST 2.2.24+\n";
	}else{
	    #program found and the correct version is installed
	}
}

=head2 fasta2stockholm

Convert a bunch of fasta files to stockholm format

=cut

sub fasta2stockholm {
	my $usage = "Usage: $0 <gapped FASTA alignment file(s)>\n";
	my @arguments = shift;
	my @argv;
	while (@arguments) {
	    my $arg = shift;
	    if ($arg =~ /^-/) {
		if ($arg eq "-h") { print $usage; exit }
		else { die $usage }
	    } else {
		push @argv, $arg;
	    }
	}
	push @argv, "-" unless @argv;

	# loop through FASTA files
	foreach my $fasta (@argv) {
	# read FASTA file
	    my %seq;
	    my @name;
	    my $name;
	    open FASTA, "<$fasta" or die "Couldn't open '$fasta': $!";
	    while (<FASTA>) {
		if (/^\s*>\s*(\S+)/) {
		    $name = $1;
		    croak "Duplicate name: $name" if defined $seq{$name};
		    push @name, $name;
		} else {
		    if (/\S/ && !defined $name) {
			warn "Ignoring: $_";
		    } else {
			s/\s//g;
			$seq{$name} .= $_;
		    }
		}
	    }
	    close FASTA;

	# check all seqs are same length
	    my $length;
	    my $lname;
	    foreach my $name (@name) {
		my $l = length $seq{$name};
		if (defined $length) {
		    croak "Sequences not all same length ($lname is $length, $name is $l)" unless $length == $l;
		} else {
		    $length = length $seq{$name};
		    $lname = $name;
		}
	    }

	# print Stockholm output
	    print "# STOCKHOLM 1.0\n";
	    foreach my $name (@name) {
		print $name, " ", $seq{$name}, "\n";
	    }
	    print "//\n";
	}
}

=head2 makeNameTable

creates a file with a table of name to marker ID mappings
Requires a marker directory as an argument

=cut

sub makeNameTable {

	my $markerDir = shift;
	`grep ">" $markerDir/*.ali > /tmp/amphora.name.table`;
	`perl -p -i -e "s/.+\:\>//g" /tmp/amphora.name.table`;
	open( NT, "/tmp/amphora.name.table" );
	while( my $line = <NT> ){
		my $commonName;
		if($line =~ /\{(.+)\}.+\[.+\]/){
			$commonName = $1;
	#		$commonName =~ s/[\.\,\/\\\(\)\:\;\'\"\{\}\$\%\^\&\*\+\-\=\s]/_/g;
		}elsif($line =~ /\[(.+)\]$/){
			$commonName = $1;
	#		$commonName =~ s/[\.\,\/\\\(\)\:\;\'\"\{\}\$\%\^\&\*\+\-\=\s]/_/g;
		}
		my @fields = split(/\s+/, $line);
		print $fields[0]."\t".$commonName."\n";
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

    perldoc Amphora2::Utilities


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

1; # End of Amphora2::Utilities
