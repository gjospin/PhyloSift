#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
BEGIN { push(@INC, "$Bin/../lib") }

use Phylosift::UpdateDB;
use Phylosift::MarkerBuild;

my @taxon_ids;
my $out_file = $ARGV[0];
my $directory = $ARGV[1];
my $marker = $ARGV[2];
my $updated = $ARGV[3];
my $pruned = $ARGV[4];
my $dna = $ARGV[5];

$Phylosift::Settings::marker_dir = $directory;
$Phylosift::Settings::ncbi_dir = "$directory/../ncbi/";
chdir $directory;
my $taxa = Phylosift::UpdateDB::filter_marker_gene_ids(marker=>$marker, updated=>$updated, pruned=>$pruned, dna=>$dna);

Phylosift::MarkerBuild::make_ncbi_subtree(out_file=>$out_file, taxon_ids=>$taxa);


