#!/usr/bin/env perl

=head1 BuildRefPack.pl

Takes in an alignment file in FASTA format as an argument and creates a reference package to use as a marker by Amphora
Usage perl BuildRefPack.pl <alignment_file> <package_name>

Amphora2/lib needs to be in the user's $PERL5PATH
WARNING : requires FastTree to be accessible from the user's $PATH

=cut

=head1 VERSION
Version 0.01

=cut

use warnings;
use strict;
use Amphora2::Utilities;
use File::Basename;

#use Amphora2::Amphora2;
use Carp;
use Cwd;
my $aln_file     = shift;
my $package_name = shift;
my $cutoff       = shift;
my $target_dir   = getcwd() . "/$package_name";
`mkdir $target_dir` unless -e $target_dir;
my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );

#my $Object = Amphora2::Amphora2->new($aln_file);
#$Object->readAmphora2Config();
#Amphora2::Utilities::programChecks();
my $hmm_file = Amphora2::Utilities::generate_hmm( $aln_file, $target_dir );

#may need to create an unaligned file for the sequences before aligning them
my $new_alignment_file = Amphora2::Utilities::hmmalign_to_model( $hmm_file, $aln_file, $target_dir );
my ( $fasttree_file, $tree_log_file ) = Amphora2::Utilities::generate_fasttree( $new_alignment_file, $target_dir );

#need to generate representatives using PDA
my $rep_file = Amphora2::Utilities::get_representatives_from_tree( $fasttree_file, $target_dir, $cutoff );

#need to read the representatives picked by PDA and generate a representative fasta file
my $rep_fasta = Amphora2::Utilities::get_fasta_from_pda_representatives( $rep_file, $target_dir, $aln_file );

#use taxit to create a new reference package required for running Amphora-2
#needed are : 1 alignment file, 1 representatives fasta file, 1 hmm profile, 1 tree file, 1 log tree file.
`cd $target_dir;taxit create -c -d "Creating a reference package for Amphora-2 for the $core marker" -l $core -f $target_dir/$core.aln -t $target_dir/$core.tree -s $target_dir/$core.log -Y FastTree -P $core`;
`rm $target_dir/$core.pda`;
`rm $target_dir/$core.tree`;
`rm $target_dir/$core.log`;
`rm $target_dir/$core.aln`;
