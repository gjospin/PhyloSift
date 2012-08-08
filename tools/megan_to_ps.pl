#! /usr/bin/perl
# Author: Eric Lowe
# Script to convert Megan output to resemble Phylosift output
# Takes a .csv file from MEGAN with the first field as read ids
# and the second field as taxids.  Must be comma separated!
# Outputs a file resembling PS output.

use strict; use warnings;

my $input = shift;
open my $fh, $input or die "Couldn't open file: $!";

while(<$fh>)
{
    chomp(my $line = $_);
    $line =~ s/ //g; # remove spaces
    $line =~ s/,/\t/g; # remove commas and replace with tabs
    print "\n$line\tno rank\tReadsreadsreads\t1\tconcat"; # print line with phylosift format
}

exit;
