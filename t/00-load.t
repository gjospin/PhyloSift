#!perl 

use Test::More tests => 1;

BEGIN {
    use_ok( 'Phylosift::Phylosift' ) || print "Bail out!\n";
}

diag( "Testing Phylosift::Phylosift $Phylosift::Phylosift::VERSION, Perl $], $^X" );
