#! /usr/bin/perl
#$ -S /usr/bin/perl 
#$ -cwd
#$ -V
# Author: Eric Lowe
# Usage: perl megan_db.pl [input file] [knockout file] [gi file]
# Script to create a custom db for MEGAN
# Extreme work in progress

use strict; use warnings;

my $kofile = shift;
my $gifile = shift;

# get an array with ko id's from subroutine
my %taxid = ();
get_koed($kofile, \%taxid);

my @dirs = qw( ebi ncbi_draft );
my $n = @dirs;

list_files(\%taxid, \@dirs, $n);
my @fnames = values %taxid;
my %koedfiles = ();
hash_kos(\@fnames, \%koedfiles);

last_push(\%koedfiles);

my %nhash = ();
hax0rz_fix(\%nhash);

get_gi($gifile, \%nhash);
write_final_hax0rz(\%nhash);

exit;


###################### SUBROUTINES ########################

# get_koed($kofile)
# sub get_koed takes in a file containing a list of taxon 
# id's that you would like to exclude from your BLAST run
# Returns an array with each element containing a separate 
# taxid.

sub get_koed
{
    my $in = shift;
    my $ref = shift;
    open my $fh, $in or die "Couldn't open knockout file $in: $!";
	
    while (<$fh>)
    {
	chomp(my $id = $_);
	if  (! defined $$ref{$id})
	{
	    $$ref{$id} = 0;
#	    print "Still good\n";
	}
    }
   
}

    
sub list_files
{
    my $ref = shift;
    my $dref = shift;
    my $n = shift;
    my $path = "/share/eisen-d2/amphora2/";
    my $dir = $$dref[$n - 1];
    $dir = $path.$dir;
#    print "$dir\n";
    my $i = 0;

    if ($n > 1)
    {
	my @files = `find $dir -name "*.fasta"`;
	foreach my $file (@files)
	{
	    chomp($file);
	    $file =~ s/$dir\///;
#	    print "$file\n";
	    $file =~ m/.*\.(\d+).fasta/;
	    my $id = $1;
#	    print "$id\n";

	    if (exists $$ref{$id})
	    {
#		print "We got one!\n";
		$$ref{$id} = $file;
		$i++;
	    }
	}
	print "$i\n";
	pop @$dref;
	$n--;
	list_files($ref, $dref,$n); 
    }elsif ($n == 1)
    {
	my @files = `find $dir -name "*.fasta"`;
        foreach my $file (@files)
        {
	    chomp($file);
            $file =~ s/$dir\///;
#	    print "$file\n";
            $file =~ m/.*\.(\d+).fasta/;
            my $id = $1;
#            print "$id\n";

            if (exists $$ref{$id})
	    {
#                print "We got one!\n";
                $$ref{$id} = $file;
                $i++;
            }
        }
        print "$i\n";
    }
}


sub hash_kos
{
    my $array = shift;
    my $ref = shift;
    my $i = 0;
    
    foreach my $file (@$array)
    {
	if (! defined $$ref{$file})
	{
	    $$ref{$file} = 1;
	    $i++;
	}
    }
    print "\$i: $i\n";

}

sub last_push
{
    my $ref = shift;
    my @dirs = qw( /share/eisen-d2/amphora2/ebi /share/eisen-d2/amphora2/ncbi_draft );
    
    foreach my $dir (@dirs)
    {
	
	my @files = `find $dir -name "*.fasta"`;
        foreach my $file (@files)
        {
	    chomp($file);
	    write_out($file, $ref);
	}
    }

}

sub write_out
{
    my $input = shift;
    my $ref = shift;

    if (exists $$ref{$input})
    {
#	print "This file got knocked the fudge out!\n";
	
    }else{
	open my $ifh, "<", $input or die "Couldn't open file for reading: $!";
	open my $ofh, ">>", "hax0rz.fasta" or die "Couldn't open your stupid hacked fasta: $!";
	my @record = ();
	while (my $line = <$ifh>)
	{
	   if ($line =~ /^>/)
	   {
	    print($ofh, @record) if @record > 1;
	    @record = ();
	   }
	   push (@record, $line);
	   chomp(my $tmp = $line);
	   next if length($tmp) < 1;
	   next unless $line =~ /\S/;
	}
	close $ifh;
	close $ofh;
    }
}

# sub hax0rz_fix
sub hax0rz_fix
{
    my $ref = shift;

    open my $fh, "<", "hax0rz.fasta" or die "Couldn't open your stupid hacked fasta: $!";

    while (<$fh>)
    {
        my $line = $_;

        if ($line =~ /^>\w+/)
        {
            if ($line =~/^>([a-zA-Z]{1,3}\d+)/)
            {
                print "Header is $1\n";
                $$ref{$1} = 1;
            }
        }
    }

    close $fh;
}

# sub get_gi
sub get_gi
{
    my $gi = shift;
    my $ref = shift;

    open my $fh, $gi or die "Couldn't open GenBank ID file $gi: $!";

    while (<$fh>)
    {
        my $line = $_;

        if ($line =~ /(\w+\d),\d,(\d+)/)
        {
            if (defined $$ref{$1})
            {
                print "Defining!\n";
                $$ref{$1} = $2;
                print "$$ref{$1} is the gi\n";
            }
        }

    }

    close $fh;
}

# sub write_final_hax0rz
sub write_final_hax0rz
{
    my $gref = shift;

    open my $fh, "<", "hax0rz.fasta" or die "Couldn't open your stupid hacked fasta: $!";
    open my $ofh, ">", "hax0rz_final.fasta" or die "Couldn't write that stupid fasta: $!";

    while (<$fh>)
    {
        my $line = $_;

        if ($line =~ /^>\w+/)
        {
            if ($line =~/^>([a-zA-Z]{1,3}\d+)/)
            {
                print "Header is $1\n";

                print "Header should be >gi\|$$gref{$1}\n";
                print $ofh ">gi|$$gref{$1}\n";
            }else{
                print $ofh "$line";
            }
        }else{
            print $ofh "$line";
        }
    }

    close $fh;
    close $ofh;
}

