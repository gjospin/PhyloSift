#!/usr/bin/perl

###Usage: pick_rep_by_mcl.pl -i input_link -n number -o output_file  -c similarity_cutoff
### if -c , -n will be invalide

use strict;
use Cwd;
my $workingDir = getcwd;
my %opt=@ARGV;
my $link_file=$opt{'-i'};
my $number_cutoff=$opt{'-n'};
if(!$number_cutoff){$number_cutoff=200;}
my $output=$opt{'-o'};

my $mcl_redunt_reduce="$workingDir/mcl_redunt_reduce.pl";
my $get_link_by_list="$workingDir/get_link_by_list.pl";

##-link test.link -include wawa -output dudu -cutoff 35 

##step1: pick a list at 95% similarity cutoff
my $cutoff=95;
my $previous_cutoff=95;
my $number_rep=0;
my $previous_number_rep=0;
my $tmp_file=$output.".reptmp1";
my $previous_tmp_file=$tmp_file.".previous";

my $log=$output.".log";
open(LOG,">$log");
print LOG "number_of_rep\tsimilarity_cutoff\n";

if(length($opt{'-c'})){$cutoff=$opt{'-c'};}

my $cmd=$mcl_redunt_reduce." -link $link_file -output $tmp_file -cutoff $cutoff";
print $cmd."\n";
system $cmd;

##if the 80 list is less than number cutoff, stop there

open(IN,$tmp_file)||die "cannot open file $tmp_file\n";
while(<IN>){
    $number_rep++;
}
close IN;


if(length($opt{'-c'}) || $number_rep <= $number_cutoff){
    rename($tmp_file,$output);
print LOG "Final: ".$number_rep."\t".$cutoff."\n";
    exit;
}

print LOG "########## ".$number_rep."\t".$cutoff."\n";

##reduce similarity cutoff at 10% interval, until the represwntative number 
## is less than the number cutoff


while($cutoff){
    $previous_cutoff=$cutoff;
    $cutoff-=10;
    rename($tmp_file,$previous_tmp_file);

$cmd=$mcl_redunt_reduce." -link $link_file -output $tmp_file -cutoff $cutoff";
print $cmd."\n";
system $cmd;

##if the 80 list is less than number cutoff, stop there
    $previous_number_rep=$number_rep;
    $number_rep=0;
open(IN,$tmp_file)||die "cannot open file $tmp_file\n";
while(<IN>){
    $number_rep++;
}
close IN;

print LOG "########## ".$number_rep."\t".$cutoff."\n";

    if($number_rep < $number_cutoff){
	last;
    }

}



$previous_number_rep=$number_rep;
$previous_cutoff=$cutoff;
rename($tmp_file,$previous_tmp_file);

## increase cutoff by 5%,

$cutoff+=5;

$cmd=$mcl_redunt_reduce." -link $link_file -output $tmp_file -cutoff $cutoff";
print $cmd."\n";
system $cmd;

    $number_rep=0;
open(IN,$tmp_file)||die "cannot open file $tmp_file\n";
while(<IN>){
    $number_rep++;
}
close IN;

print LOG "########## ".$number_rep."\t".$cutoff."\n";

## what I get is less than the number cutoff, increase the cutoff by 1% until it is more than the number cutoff
## if what I get is more than the number cutoff, the bottom is the previous cutoff

my $bottom_cutoff;
if($number_rep < $number_cutoff){
$bottom_cutoff=$cutoff;
rename($tmp_file,$previous_tmp_file);
}
else{
$bottom_cutoff=$previous_cutoff;
unlink $tmp_file;
}
 


for ($cutoff=$bottom_cutoff+1;$cutoff<=$bottom_cutoff+3;$cutoff++){

$cmd=$mcl_redunt_reduce." -link $link_file -output $tmp_file -cutoff $cutoff";
print $cmd."\n";
system $cmd;

$number_rep=0;
open(IN,$tmp_file)||die "cannot open file $tmp_file\n";
while(<IN>){
    $number_rep++;
}
close IN;

print LOG "########## ".$number_rep."\t".$cutoff."\n";

if($number_rep > $number_cutoff){last;}

rename ($tmp_file,$previous_tmp_file);
$previous_cutoff=$cutoff;
$previous_number_rep=$number_rep;
}


rename ($previous_tmp_file,$output);
unlink $tmp_file;


print LOG "Final ".$previous_number_rep."\t".$previous_cutoff."\n";
close LOG;
