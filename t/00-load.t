#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Amphora2::Amphora2' ) || print "Bail out!\n";
}

diag( "Testing Amphora2::Amphora2 $Amphora2::Amphora2::VERSION, Perl $], $^X" );
