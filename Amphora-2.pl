#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Amphora2::Amphora2;
my $pair=0;
my $threadNum=1;
my $clean=0;
my $force=0;
my $continue=0;
my $isolate=0;
my $custom = "";
my $usage = qq~
  Usage: $0 <mode> <options> <reads_file>                                                                       

    ~;
my $usage2 = qq~
  Usage: $0 <mode> <options> -paired <reads_file_1> <reads_file_2>                                              

    ~;


GetOptions("threaded=i" => \$threadNum,
	   "clean" => \$clean,
	   "paired" => \$pair, # used for paired fastQ input split in 2 different files                           
	   "custom=s" => \$custom, #need a file containing the marker names to use without extensions ** marker names shouldn't contain '_'
	   "f" => \$force, #overrides a previous run otherwise stop                                               
	   "continue" => \$continue, #when a mode different than all is used, continue the rest of Amphora after the section
	   "isolate" => \$isolate, #use when processing one or more isolate genomes
	   #specified by the mode is finished
    )|| die $usage;
    

print "@ARGV\n";
if($pair ==0){
    croak $usage unless ($ARGV[0] && $ARGV[1]);
}elsif(scalar(@ARGV) == 3){
    croak $usage2 unless ($ARGV[0] && $ARGV[1] && $ARGV[2]);
}else{
    die $usage."\nOR\n".$usage2."\n";
}



my $newObject = new Amphora2::Amphora2();
print "PAIR : $pair\n";
if($pair == 0){
    print "notpaired @ARGV\n";
    $newObject = $newObject->initialize($ARGV[0],$ARGV[1]);
}else{
    print "INSIDE paired\n";
    $newObject = $newObject->initialize($ARGV[0],$ARGV[1],$ARGV[2]);
}

$newObject->{"isolate"} = $isolate;

my $readsFile = $newObject->getReadsFile;
print "FORCE: ".$force."\n";
print "Continue : ".$continue."\n";
$newObject->run($force,$custom,$continue,$isolate);

#Amphora2::Amphora2::run(@ARGV);

