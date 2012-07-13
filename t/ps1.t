#!perl
# Test module for Phylosift::Phylosift
# Need to add class/object testing and Test::Harness
# Look at adding Devel::Cover as well for 'busy parts' of code

use strict; use warnings;

use Test::More 'no_plan';

my @subs = qw( initialize getReadsFile run read_phylosift_config run_program_check );

use_ok( 'Phylosift::Phylosift', @subs );

use_ok( 'Phylosift::Phylosift' );

my $ps = Phylosift->new( );
isa_ok( $ps, 'Phylosift::Phylosift' );
