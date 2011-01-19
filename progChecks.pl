#!/usr/bin/perl -w

# this script checks the program requirements for Amphora-2
# writes to STDERR if a program is missing or the wrong version is installed
# returns 1 or 0 depending on success of failure.

use strict;
use warnings;


eval 'require Bio::Seq;';
if ($@) {
    print STDERR "Bioperl was NOT found\n";
}


my $pplacercheck = `which pplacer`;
if($pplacercheck eq ""){
    #program not found return;
    exit(1);
    die("pplacer v1.1.alpha00 not found");
}elsif(`pplacer -help` !~ m/v1.1.alpha00/){
    # pplacer was found but the version doens't match the one tested with Amphora
    print STDERR "Warning : a different version of pplacer was found. Amphora-2 was tested with pplacer v1.1.alpha00\n";
}else{
    #program found and the correct version is installed
}

my $hmmercheck = `which hmmalign`;
if($hmmercheck eq ""){
    #program not found return;
    exit(1);
    die("HMMER3 not found");
}elsif(`hmmalign -h` !~ m/HMMER 3.0rc1/){
    # pplacer was found but the version doens't match the one tested with Amphora
    print STDERR "Warning : a different version of HMMER was found. Amphora-2 was tested with HMMER 3.0rc1\n";
}else{
    #program found and the correct version is installed
}


my $blastcheck = `which blastp`;
if($blastcheck eq ""){
    #program not found return;
    exit(1);
    die("Blast 2.2.24+ not found");
}elsif(`blastp -help` !~ m/BLAST 2.2.24+/){
    # pplacer was found but the version doens't match the one tested with Amphora
    print STDERR "Warning : a different version of Blast was found. Amphora-2 was tested with BLAST 2.2.24+\n";
}else{
    #program found and the correct version is installed
}
exit(0);
