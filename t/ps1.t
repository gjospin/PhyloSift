#!perl

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

ok( placeholder_1() == 4 );
ok( placeholder_2() eq "Working?" );
