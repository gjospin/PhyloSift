#!perl
# Test module for Phylosift::Phylosift
# Need to add class/object testing and Test::Harness
# Look at adding Devel::Cover as well for 'busy parts' of code

use strict; use warnings;

use Test::More 'no_plan';

sub placeholder_1
{
	my $i = 2 + 2;
	return $i;
}

sub placeholder_2
{
    my $check = "Working?";
    return $check;
}

my @subs = qw( initialize getReadsFile run read_phylosift_config run_program_check );

use_ok( 'Phylosift::Phylosift', @subs );

ok( placeholder_1() == 4, 'placeholder_1() output should be 4' );
ok( placeholder_2() eq "Working?", 'placeholder_2() output should be sane' );
