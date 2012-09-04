#! /usr/bin/perl
# Author: Eric Lowe
use strict; 
use warnings;

use Test::More qw(no_plan);
#use Test::Warn;
use lib '../lib';

my @subs = qw(ps_open);
BEGIN { use_ok( 'Phylosift::Utilities', @subs ) or exit; }
require_ok('Phylosift::Utilities') or exit;

my $file = "testfile.txt";

is(ps_open("<", $file), $file, "ps_open");

TODO: {
    local $TODO = 'Working on testing warnings';
    warnings_are{is(ps_open("<", $file), $file, "ps_open")} ['Unable to read from testfile.txt'];
}
# done_testing();