#!/usr/bin/env perl
#
# script to download and package up a standalone version of phylosift
#
use strict;
use warnings;
`rm -rf PhyloSift`;

my $branch = $ARGV[0] || "master";

#`rm -rf bioperl-live`;
`git clone git://github.com/gjospin/PhyloSift.git`;
`cd PhyloSift ; git checkout $branch`;
`rm -rf PhyloSift/.git`;
`rm PhyloSift/Makefile.PL`;
`rm PhyloSift/ignore.txt`;
`rm PhyloSift/MANIFEST`;
`rm PhyloSift/Changes`;
`rm PhyloSift/.gitignore`;
`rm PhyloSift/.includepath`;
`rm PhyloSift/.project`;
`rm -rf PhyloSift/t`;
`rm -rf PhyloSift/tools`;

# add bioperl-live
`git clone git://github.com/bioperl/bioperl-live.git`;
`mv bioperl-live/Bio* PhyloSift/lib`;

# add Rutger Vos' Bio::Phylo
`curl -LO http://search.cpan.org/CPAN/authors/id/R/RV/RVOSA/Bio-Phylo-0.45.tar.gz`;
`tar xvzf Bio-Phylo-0.45.tar.gz`;
chdir("Bio-Phylo-0.45");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Bio/Phylo* ../PhyloSift/lib/Bio/`;
chdir("..");

# add JSON package
`curl -LO http://search.cpan.org/CPAN/authors/id/M/MA/MAKAMAKA/JSON-2.53.tar.gz`;
`tar xvzf JSON-2.53.tar.gz`;
chdir("JSON-2.53");
`perl Makefile.PL`;
`make`;
`mv blib/lib/JSON* ../PhyloSift/lib/`;
chdir("..");

# Encode::Locale
`curl -LO http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/Encode-Locale-0.04.tar.gz`;
`tar xzf Encode-Locale-0.04.tar.gz`;
chdir("Encode-Locale-0.04");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Encode ../PhyloSift/lib/`;
chdir("..");

# add Locale::Maketext
`curl -LO http://search.cpan.org/CPAN/authors/id/T/TO/TODDR/Locale-Maketext-1.19.tar.gz`;
`tar xvzf Locale-Maketext-1.19.tar.gz`;
chdir("Locale-Maketext-1.19");

# remove the following files because they break Todd's ancient perldoc
`rm lib/Locale/Maketext/*.pod`;
`perl Makefile.PL`;
`make`;
`mv blib/lib/Locale/ ../PhyloSift/lib/`;
chdir("..");

# Locale::Maketext::Simple 
`curl -LO http://search.cpan.org/CPAN/authors/id/J/JE/JESSE/Locale-Maketext-Simple-0.21.tar.gz`;
`tar xvzf Locale-Maketext-Simple-0.21.tar.gz`;
chdir("Locale-Maketext-Simple-0.21");
`perl Makefile.PL`;
`make`;
`mv blib/lib/Locale/Maketext/* ../PhyloSift/lib/Locale/Maketext/`;
chdir("..");

# XML::Writer
`curl -LO http://search.cpan.org/CPAN/authors/id/J/JO/JOSEPHW/XML-Writer-0.615.tar.gz`;
`tar xzf XML-Writer-0.615.tar.gz`;
chdir("XML-Writer-0.615");
`perl Makefile.PL`;
`make`;
`mv blib/lib/XML/ ../PhyloSift/lib/`;
chdir("..");

# libwww-perl (for LWP::Simple)
`curl -LO http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/libwww-perl-6.04.tar.gz`;
`tar xzf libwww-perl-6.04.tar.gz`;
chdir("libwww-perl-6.04");
`perl Makefile.PL`;
`make`;
`mv blib/lib/LWP/ ../PhyloSift/lib/`;
chdir("..");

`curl -LO http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/HTTP-Message-6.03.tar.gz`;
`tar xzf HTTP-Message-6.03.tar.gz`;
chdir("HTTP-Message-6.03");
`perl Makefile.PL`;
`make`;
`mv blib/lib/HTTP/ ../PhyloSift/lib/`;
chdir("..");

`curl -LO http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/HTTP-Date-6.02.tar.gz`;
`tar xzf HTTP-Date-6.02.tar.gz`;
chdir("HTTP-Date-6.02");
`perl Makefile.PL`;
`make`;
`mv blib/lib/HTTP/ ../PhyloSift/lib/`;
chdir("..");


# add Version.pm
`curl -LO http://search.cpan.org/CPAN/authors/id/J/JP/JPEACOCK/version-0.95.tar.gz`;
`tar xzf version-0.95.tar.gz`;
chdir("version-0.95");
`perl Makefile.PL --perl_only`;
`make`;

# put these in "legacy" because we only want to use them if the perl version is ancient -- including them breaks newer perls
`mv blib/lib/version* ../PhyloSift/legacy/`;

chdir("..");

# package everything up and datestamp it
my @timerval = localtime();
my $datestr  = ( 1900 + $timerval[5] );
$datestr .= 0 if $timerval[4] <= 9;
$datestr .= ( $timerval[4] + 1 );
$datestr .= 0 if $timerval[3] <= 9;
$datestr .= $timerval[3];
`mv PhyloSift phylosift_$datestr`;
`tar cjf phylosift_$datestr.tar.bz2 phylosift_$datestr`;
`rm -rf phylosift_$datestr`;
`echo "phylosift_$datestr" > psversion`;
exit 0;
