#! /usr/bin/perl
# Author: Eric Lowe
use strict; 
use warnings;

use Test::More qw(no_plan);

BEGIN { use_ok( 'Phylosift::Utilities' ); }
require_ok('Phylosift::Utilities');

my $mode = "<";
my $file = "testfile.txt";

is(ps_open($mode, $file), $file, "ps_open");

# done_testing();