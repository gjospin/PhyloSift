#!/usr/bin/env perl
#
# script to download and package up a standalone version of amphora2
#
use strict;
use warnings;
`rm -rf Amphora-2`;
#`rm -rf bioperl-live`;
`git clone git://github.com/gjospin/Amphora-2.git`;
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
`tar czf Amphora-2.tar.gz Amphora-2`;

