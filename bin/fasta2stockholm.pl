#!/usr/bin/env perl -w

my $usage = "Usage: $0 <gapped FASTA alignment file(s)>\n";

my @argv;
while (@ARGV) {
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
	    die "Duplicate name: $name" if defined $seq{$name};
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
	    die "Sequences not all same length ($lname is $length, $name is $l)" unless $length == $l;
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

