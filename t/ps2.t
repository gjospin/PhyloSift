#!perl
# Test module for Phylosift::Utilities
# As with ps1.t, will add Test::Harness and Devel::Cover
use strict; use warnings;

use Test::More tests => 2;
use_ok( 'Phylosift::Utilities' ) or exit;
use_ok( 'Phylosift::Summarize' ) or exit;

