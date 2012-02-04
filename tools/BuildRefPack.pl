#!/usr/bin/env perl


=head1 BuildRefPack.pl

Takes in an alignment file in FASTA format as an argument and creates a reference package to use as a marker by Amphora
Usage perl BuildRefPack.pl <alignment_file> <package_name>

WARNING : requires FastTree to be accessible from the user's $PATH

=cut

=head1 VERSION

Version 0.01

=cut

use warnings;
use strict;
use Amphora2::Utilities;
#use Amphora2::Amphora2;
use Carp;
use Cwd;

my $aln_file = shift;
my $package_name = shift;
my $target_dir = getcwd()."/$package_name";
`mkdir $target_dir` unless -e $target_dir;

#my $Object = Amphora2::Amphora2->new($aln_file);
#$Object->readAmphora2Config();
#Amphora2::Utilities::programChecks();
my $hmm_file = Amphora2::Utilities::generate_hmm($aln_file,$target_dir);
#may need to create an unaligned file for the sequences before aligning them
my $new_alignment_file = Amphora2::Utilities::hmmalign_to_model($hmm_file,$aln_file,$target_dir);
my $fasttree_file = Amphora2::Utilities::generate_fasttree($new_alignment_file,$target_dir);

