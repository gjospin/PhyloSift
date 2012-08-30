#! /usr/bin/perl
# Test module for Phylosift::Phylosift
# Need to add class/object testing and Test::Harness
# Look at adding Devel::Cover as well for 'busy parts' of code

use strict;
use warnings;

use Test::More qw(no_plan);
#use Test::Warn;
BEGIN { chdir 't' if -d 't' }
use lib '../lib';
use blib;

BEGIN { use_ok( 'Phylosift::Phylosift' ) or exit; }
require_ok('Phylosift::Phylosift');

TODO: {
	local $TODO = 'Still working out some kinks';
    my $ps = Phylosift->new( );
    isa_ok( $ps, 'Phylosift::Phylosift' );
}
