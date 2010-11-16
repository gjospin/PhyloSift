#!/usr/bin/perl

##take in a link file established by /dwu_scripts/public_scripts/public_scripts/blastall2link.pl, output an nonredundant list:
## the minimun similarity cutoff
  

use strict;


my %opt=@ARGV;
my $link_file=$opt{'-link'};
my $list_file=$opt{'-include'};
my $cutoff=$opt{'-cutoff'};
my $output=$opt{'-output'};


if(!(length($cutoff)>=1)){$cutoff=0;}
my %in;


    open(IN,$list_file)||die "cannot open file $list_file\n";
    while(<IN>){
	my @t=split(/\s+/);
	foreach my $t(@t){
	    $in{$t}=1;
	}
    }
    close IN;


open(OUT,">$output") || die "cannot create file $output \n";
open(IN, $link_file)||die "cannot open file $link_file\n";
while(<IN>){
    my ($acc1,$acc2,$similarity)=split(/\s+/);
    if(!($in{$acc1}&&$in{$acc2})){next;}
    if($similarity < $cutoff){next;}
    print OUT $acc1."\t".$acc2."\t".$similarity."\n";
}
close IN;
close OUT;

