#!/usr/bin/env perl
#
# script to download and package up a standalone version of phylosift
#
use strict;
use warnings;
`rm -rf PhyloSift`;

#`rm -rf bioperl-live`;
`git clone git://github.com/gjospin/PhyloSift.git`;
`rm -rf PhyloSift/.git`;
`git clone git://github.com/bioperl/bioperl-live.git`;
`mv bioperl-live/Bio* PhyloSift/lib`;
`wget http://search.cpan.org/CPAN/authors/id/R/RV/RVOSA/Bio-Phylo-0.45.tar.gz`;
`tar xvzf Bio-Phylo-0.45.tar.gz`;
chdir("Bio-Phylo-0.45");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Bio/Phylo* ../PhyloSift/lib/Bio/`;
chdir("..");
`wget http://search.cpan.org/CPAN/authors/id/T/TO/TODDR/Locale-Maketext-1.19.tar.gz`;
`tar xvzf Locale-Maketext-1.19.tar.gz`;
chdir("Locale-Maketext-1.19");

# remove the following files because they break Todd's ancient perldoc
`rm lib/Locale/Maketext/*.pod`;
`perl Makefile.PL`;
`make`;
`mv blib/lib/Locale/ ../PhyloSift/lib/`;
chdir("..");
`wget http://search.cpan.org/CPAN/authors/id/J/JP/JPEACOCK/version-0.95.tar.gz`;
`tar xzf version-0.95.tar.gz`;
chdir("version-0.95");
`perl Makefile.PL`;
`make`;

# put these in "legacy" because we only want to use them if the perl version is ancient -- including them breaks newer perls
`mkdir -p ../PhyloSift/legacy/arch/auto`;
`mv blib/lib/version* ../PhyloSift/legacy/`;
`mv blib/arch/auto/version ../PhyloSift/legacy/arch/auto/`;
chdir("..");
my @timerval = localtime();
my $datestr  = ( 1900 + $timerval[5] );
$datestr .= 0 if $timerval[4] < 9;
$datestr .= ( $timerval[4] + 1 );
$datestr .= 0 if $timerval[3] < 9;
$datestr .= $timerval[3];
`mv PhyloSift phylosift_$datestr`;
`tar cjf phylosift_$datestr.tar.bz2 phylosift_$datestr`;
`rm -rf phylosift_$datestr`;
`cp phylosift_$datestr.tar.bz2 ~/public_html/phylosift/phylosift_latest.tar.bz2`;
`mv phylosift_$datestr.tar.bz2 ~/public_html/phylosift/`;
