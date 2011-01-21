#!/usr/bin/perl -w
use strict;
my $markerDir = $ARGV[0];
`grep ">" $markerDir*.ali > /tmp/amphora.name.table`;
`perl -p -i -e "s/.+\:\>//g" /tmp/amphora.name.table`;
open( NT, "/tmp/amphora.name.table" );
while( my $line = <NT> ){
	my $commonName;
	if($line =~ /\[(.+)\]$/){
		$commonName = $1;
		$commonName =~ s/[\.\,\/\\\(\)\:\;\'\"\{\}\$\%\^\&\*\+\-\=\s]/_/g;
	}
	my @fields = split(/\s+/, $line);
	print $fields[0]."\t".$commonName."\n";
}
