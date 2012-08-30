#! /usr/bin/perl
# Author: Eric Lowe
use strict; 
use warnings;

use Test::More qw(no_plan);
my @subs = qw(ps_open);
BEGIN { use_ok( 'Phylosift::Utilities', @subs ) or exit; }
require_ok('Phylosift::Utilities') or exit;

my $file = "testfile.txt";

is(ps_open("<", $file), $file, "ps_open");

# done_testing();