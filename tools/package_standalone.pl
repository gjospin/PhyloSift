#!/usr/bin/env perl
#
# script to download and package up a standalone version of amphora2
#
use strict;
use warnings;
`rm -rf Amphora-2`;
#`rm -rf bioperl-live`;
`git clone git://github.com/gjospin/Amphora-2.git`;
`rm -rf Amphora-2/.git`;
`git clone git://github.com/bioperl/bioperl-live.git`;
`mv bioperl-live/Bio* Amphora-2/lib`;
`wget http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.71.tar.gz`;
`tar xvzf Math-Random-0.71.tar.gz`;
chdir("Math-Random-0.71");
`perl Makefile.PL`;
`make`;
`mv blib/arch/auto ../Amphora-2/lib/`;
`mv blib/lib/Math ../Amphora-2/lib/`;
chdir("..");

`wget http://search.cpan.org/CPAN/authors/id/R/RV/RVOSA/Bio-Phylo-0.45.tar.gz`;
`tar xvzf Bio-Phylo-0.45.tar.gz`;
chdir("Bio-Phylo-0.45");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Bio/Phylo* ../Amphora-2/lib/Bio/`;
chdir("..");

`wget http://cpan.uwinnipeg.ca/cpan/authors/id/T/TO/TODDR/Locale-Maketext-1.19.tar.gz`;
`tar xvzf Locale-Maketext-1.19.tar.gz`;
chdir("Locale-Maketext-1.19");
# remove the following files because they break Todd's ancient perldoc
`rm lib/Locale/Maketext/*.pod`;
`perl Makefile.PL`;
`make`;
`mv blib/lib/Locale/ ../Amphora-2/lib/`;
chdir("..");

`wget http://search.cpan.org/CPAN/authors/id/J/JP/JPEACOCK/version-0.95.tar.gz`;
`tar xzf version-0.95.tar.gz`;
chdir("version-0.95");
`perl Makefile.PL`;
`make`;
# put these in "legacy" because we only want to use them if the perl version is ancient -- including them breaks newer perls
`mkdir -p ../Amphora-2/legacy/arch/auto`;
`mv blib/lib/version* ../Amphora-2/legacy/`;
`mv blib/arch/auto/version ../Amphora-2/legacy/arch/auto/`;
chdir("..");

my @timerval = localtime();
my $datestr = (1900+$timerval[5]);
$datestr .= 0 if $timerval[4] < 9; 
$datestr .= ($timerval[4]+1);
$datestr .= 0 if $timerval[3] < 9; 
$datestr .= $timerval[3];

`mv Amphora-2 amphora2_$datestr`;
`tar cjf amphora2_$datestr.tar.bz2 amphora2_$datestr`;
`rm -rf amphora2_$datestr`;
`cp amphora2_$datestr.tar.bz2 ~/public_html/amphora2/amphora2_latest.tar.bz2`;
`mv amphora2_$datestr.tar.bz2 ~/public_html/amphora2/`;
