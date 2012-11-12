#! /usr/bin/perl
# Author: Eric Lowe
# Settings.pm is a very small module
# This test should check its one subroutine effectively
use strict;
use warnings;

use Test::More tests => 2;
use lib '../lib';

my @subs = qw(set_default);
BEGIN { use_ok( 'Phylosift::Settings', @subs ) or exit; }

my $testParam;
my $expected = 5;
# should set default value of "parameter" to 5
set_default(parameter=> \$testParam, value=> 5);
# tests if parameter is 5 by cmp_ok test
cmp_ok($testParam, '==', $expected, 'value == 5');

