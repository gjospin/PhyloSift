#! /usr/bin/perl
#$ -S /usr/bin/perl 
#$ -cwd
#$ -V
# Author: Eric Lowe
# Usage: perl megan_db.pl [knockout file] [gi file]
# Script to create a custom db for MEGAN
# Extreme work in progress

use strict; use warnings;

my $kofile = shift; # file of knockout taxids
my $gifile = shift; # file of gb accession to gi accession

# get a hash with ko id's from subroutine
my %taxid = ();
get_koed($kofile, \%taxid);
# list of directories to search for files
my @dirs = qw( ebi ncbi_draft );
my $n = @dirs; # number of directories in list

list_files(\%taxid, \@dirs, $n); # get the list of appropriate files
my @fnames = values %taxid; # array of filenames we don't want to include
my %koedfiles = (); # hash of files to be knocked out
hash_kos(\@fnames, \%koedfiles); # defines a hash with keys of files to be knocked out

last_push(\%koedfiles); # begins the first writing out, calls a subroutine that checks if file
                        # is to be used or not, writes to hax0rz.fasta
                        
my %nhash = (); # hash of gb accessions and related gi accessions to replace headers
hax0rz_fix(\%nhash); # finds headers to fix

get_gi($gifile, \%nhash); # finds appropriate gi to replace gb accession
write_final_hax0rz(\%nhash); # writes final file with correct headers

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
	    print $ofh @record if @record > 1;
	    @record = ();
	   }
	   push (@record, $line);
	   chomp(my $tmp = $line);
	   next if length($tmp) < 1;
	   next unless $line =~ /\S/;
	}
	print $ofh @record if @record > 1; # prints anything left to output file
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
    my %checkgi = ();
    my $defined = 0;
    
    open my $fh, "<", "hax0rz.fasta" or die "Couldn't open your stupid hacked fasta: $!";
    open my $ofh, ">", "hax0rz_final.fasta" or die "Couldn't write that stupid fasta: $!";

    while (<$fh>)
    {
        my $line = $_;

        if ($line =~ /^>\w+/)
	{
            if ($line =~ /^>([a-zA-Z]{1,3}\d+)/)
            {
             	print "Header is $1\n";
                my $g = $$gref{$1};

                print "Header should be >gi\|$g\n";
		if (! defined $checkgi{$g})
                {
                    $checkgi{$g} = 1;
                    $defined = 0;
                    print $ofh ">gi|$g\n";
                }else{
                    $defined = 1;
                }
            }elsif ($line =~ /^>gi\|(\d+)/)
            {
                if (! defined $checkgi{$1})
                {
                    $checkgi{$1} = 1;
                    $defined = 0;
                    print $ofh "$line";
                }else{
                    $defined = 1;
                }
            }
        }else{
            print $ofh "$line" unless $defined == 1;
        }
    }

    close $fh;
    close $ofh;
}

