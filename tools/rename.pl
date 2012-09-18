#! /usr/bin/perl
# Author: Eric Lowe
# This is a simple Perl script to facilitate building
# on Buildbot buildslaves
# fixes issue with buildslaves not unpacking phylosift tarball

use strict; use warnings;

open FILE, "< psversion" or die "Couldn't find psversion!";
chomp(my $name = <FILE>);      # Get Phylosift tarball name from psversion
my $tarname = $name. '.tar.bz2';
`tar -xjf $tarname`;
`mv $name ../phylosift`;


exit;
