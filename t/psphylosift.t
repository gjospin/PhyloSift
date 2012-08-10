#! /usr/bin/perl
# Test module for Phylosift::Phylosift
# Need to add class/object testing and Test::Harness
# Look at adding Devel::Cover as well for 'busy parts' of code

use strict;
use warnings;

use Test::More qw(no_plan);

BEGIN { use_ok( 'Phylosift::Phylosift' ); }
require_ok('Phylosift::Phylosift');

#my $ps = Phylosift->new( );
#isa_ok( $ps, 'Phylosift::Phylosift' );