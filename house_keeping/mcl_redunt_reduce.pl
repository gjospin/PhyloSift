#!/usr/bin/perl

##take in a link file established by /dwu_scripts/public_scripts/public_scripts/blastall2link.pl, output an nonredundant list:
## the minimun similarity cutoff
  
##requires : 
##   /home/dywu/bin/bin/mcl
## Usage: mcl_redunt_reduce.pl -link link_file -include list_of_accs 
##        -cutoff cutoff -output output_file

use strict;
my $mcl="/home/gjospin/programs/bin/mcl";
my %opt=@ARGV;
my $link_file=$opt{'-link'};
my $list_file=$opt{'-include'};
my $cutoff=$opt{'-cutoff'};
my $output=$opt{'-output'};
my $tmp_link=$output.".tmplink";
my $tmp_cluster=$output.".tmpcluster";

if(!(length($cutoff)>=1)){$cutoff=0;}
my $part=0;
my %in;

if(-s $list_file){
    $part=1;
    open(IN,$list_file)||die "cannot open file $list_file\n";
    while(<IN>){
	my @t=split(/\s+/);
	foreach my $t(@t){
	    $in{$t}=1;
	}
    }
    close IN;
}


##print "cutoff $cutoff \n";
open(OUT,">$tmp_link") || die "cannot create file $tmp_link \n";
open(IN, $link_file)||die "cannot open file $link_file\n";
while(<IN>){
    my ($acc1,$acc2,$similarity)=split(/\s+/);
    if($part && !($in{$acc1}&&$in{$acc2})){next;}
    if($similarity < $cutoff){next;}
    print OUT $acc1."\t".$acc2."\t".$similarity."\n";
}
close IN;
close OUT;

my $cmd=$mcl." $tmp_link --abc -I 2.0 -o $tmp_cluster";

system $cmd;

open(IN,$tmp_cluster)||die "cannot open file $tmp_cluster\n";
open(OUT,">$output") || die "cannot create file $output\n";
while(<IN>){
    my ($acc)=split(/\s+/);
    print OUT $acc."\n";
}

close OUT;
close IN;

unlink $tmp_link;
unlink $tmp_cluster;
