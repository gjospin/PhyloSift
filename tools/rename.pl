#! /usr/bin/perl
# Author: Eric Lowe
# This is a simple Perl script to facilitate building
# on Buildbot buildslaves
# fixes issue with buildslaves not unpacking phylosift tarball

use strict; use warnings;

open FILE, "< psversion" or die "Couldn't find psversion!";
chomp(my $name = <FILE>);      # Get Phylosift tarball name from psversion
$name = $name . '.tar.bz2';    # Append extension
my $newname = 'phylosift.tar.bz2';   # new name for tarball
rename ($name, $newname);      # rename files

exit;
