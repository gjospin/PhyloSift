#!/usr/bin/perl -w
use strict;

if(@ARGV!=1){
	print STDERR "Usage: summarize.pl <amphora output>\n";
}

my $markerdir = "markers";
open( NAMETABLE, "$markerdir/name.table" );
my %namemap;
while( my $line = <NAMETABLE> ){
	chomp $line;
	my @vals = split( /\t/, $line );
	$namemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]);
}

open( TAXIDS, "ncbi/names.dmp" );
my %nameidmap;
my %idnamemap;
while( my $line = <TAXIDS> ){
	chomp $line;
	if(($line =~ /scientific name/) || ($line =~ /synonym/) || ($line =~ /misspelling/)){
		my @vals = split( /\s+\|\s+/, $line );
		$nameidmap{homogenizeNameAlaDongying($vals[1])}=$vals[0];
		$idnamemap{$vals[0]}=homogenizeNameAlaDongying($vals[1]) if($line =~ /scientific name/);
	}
}

open( TAXSTRUCTURE, "ncbi/nodes.dmp" );
my @parent;
while( my $line = <TAXSTRUCTURE> ){
	chomp $line;
	my @vals = split( /\s+\|\s+/, $line );
	$parent[$vals[0]] = $vals[1];
}

open( TAXACALLS, $ARGV[0] );

my %hitcounter;
my $readcount = 0;
while( my $line = <TAXACALLS> ){
	chomp $line;
	$line =~ s/\s+$//g;
	$line =~ s/^\s+//g;
#	next unless length($line) > 12;
	next unless( defined($namemap{$line}) );
	my ($tid,$name) = dongyingFindNameInTaxaDb($namemap{$line});
	if($tid eq "ERROR"){
		print STDERR "Error! Could not find $line in name map\n" if length($line) > 12;
		next;
	}
	
	#got the taxon id, now walk to root tallying everything we hit
	next unless(defined($tid));
	while( $tid != 1 ){
		if(defined($hitcounter{$tid})){
			$hitcounter{$tid}++;
		}else{
			$hitcounter{$tid}=1;
		}
		$tid = $parent[$tid];
	}
	$readcount++;
}

my %hitvals;
foreach my $tid(keys(%hitcounter)){
	my $frac = $hitcounter{$tid}/$readcount;
	$hitvals{$idnamemap{$tid}} = $frac;
}
my @sorted = reverse sort { $hitvals{$a} cmp $hitvals{$b} } keys %hitvals; 

foreach my $names(@sorted){
	print "$names\t".$hitvals{$names}."\n";
}


sub homogenizeNameAlaDongying {
	my $inName = shift;
	$inName=~s/^\s+//;
	$inName=~s/\s+$//;
	$inName=~s/\s+/ /g;
	$inName=~s/,//g;
	$inName=uc $inName;
	return $inName;
}

sub dongyingFindNameInTaxaDb {
	my $name = shift;
	$name=~s/^\s+//;
	my @t=split(/\s+/, $name);
	my $input_name=join(" ",@t);
	my $q_name=$input_name;
	my $id="ERROR";
	while(@t>=1){
		$q_name=join(" ",@t);
		$q_name=uc $q_name;
		if(defined($nameidmap{$q_name})){
			$id=$nameidmap{$q_name};
			last;
		}
		pop @t;
	}
	return ($id,$q_name);
}

