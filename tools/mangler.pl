#!/usr/bin/env perl
use strict;
use warnings;

while( my $line = <STDIN> ){
	$line =~ s/:(.+?)\{(\d+?)\}/\{$2\}:$1/g;
	print $line;
}

