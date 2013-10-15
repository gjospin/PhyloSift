#!/usr/bin/env perl
#
# script to download and package up a standalone version of phylosift
#
use strict;
use warnings;
use File::Basename;

my $prefix = "../PhyloSift";

if (@ARGV>0 && $ARGV[0] eq "buildbot") # if buildbot
{
	$prefix = ".."; # set path prefix
	
	# add bioperl latest release (working in buildbot)
    `git clone git://github.com/bioperl/bioperl-live.git`;
    `cd bioperl-live ; git checkout 1d6bccb6d121f050de44dd86b80a02d6abb159a3; cd ..`;
    `mv bioperl-live/Bio* lib/`;
}else{
	`rm -rf PhyloSift`;
	
	my $branch = $ARGV[0] || "master"; # set branch to checkout
    $prefix = "../PhyloSift"; # set path prefix
    
    `rm -rf bioperl-live`;
    `git clone https://github.com/gjospin/PhyloSift.git`;
    `cd PhyloSift ; git checkout $branch`;
    `rm -rf PhyloSift/.git`;
    `rm PhyloSift/Makefile.PL`;
    `rm PhyloSift/ignore.txt`;
    `rm PhyloSift/MANIFEST`;
    `rm PhyloSift/Changes`;
    `rm PhyloSift/.gitignore`;
    `rm PhyloSift/.includepath`;
    `rm PhyloSift/.project`;
    `rm -rf PhyloSift/t`;
    `rm -rf PhyloSift/tools`;
    `rm -rf PhyloSift/web`;
    
    # add bioperl latest release (not on buildbot)
    `git clone https://github.com/bioperl/bioperl-live.git`;
    `cd bioperl-live ; git checkout 1d6bccb6d121f050de44dd86b80a02d6abb159a3; cd ..`;
    `mv bioperl-live/Bio* PhyloSift/lib/`;
}

# add Rutger Vos' Bio::Phylo
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RV/RVOSA/Bio-Phylo-0.56.tar.gz", mv_cmd=>"mv lib/Bio/* $prefix/lib/Bio/");

# add JSON package
add_package(url=>"http://search.cpan.org/CPAN/authors/id/M/MA/MAKAMAKA/JSON-2.53.tar.gz");

# Encode::Locale
add_package(url=>"http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/Encode-Locale-1.03.tar.gz");

# add Locale::Maketext
`curl -LO http://search.cpan.org/CPAN/authors/id/T/TO/TODDR/Locale-Maketext-1.19.tar.gz`;
`tar xvzf Locale-Maketext-1.19.tar.gz`;
chdir("Locale-Maketext-1.19");

# remove the following files because they break Todd's ancient perldoc
`rm lib/Locale/Maketext/*.pod`;
`perl Makefile.PL`;
`make`;
`mv blib/lib/Locale/ $prefix/lib/`;
chdir("..");

# Locale::Maketext::Simple 
add_package(url=>"http://search.cpan.org/CPAN/authors/id/J/JE/JESSE/Locale-Maketext-Simple-0.21.tar.gz", mv_cmd=>"mv blib/lib/Locale/Maketext/* $prefix/lib/Locale/Maketext/");

# XML::Writer
add_package(url=>"http://search.cpan.org/CPAN/authors/id/J/JO/JOSEPHW/XML-Writer-0.615.tar.gz");

# File::NSFLock
add_package(url=>"http://search.cpan.org/CPAN/authors/id/B/BB/BBB/File-NFSLock-1.21.tar.gz");

# add Version.pm
#`mkdir -p $prefix/legacy`; #prefix is not right here as we are in the top level directory
`mkdir -p PhyloSift/legacy`;
add_package(url=>"http://search.cpan.org/CPAN/authors/id/J/JP/JPEACOCK/version-0.9902.tar.gz", make_opts=>"--perl_only", mv_cmd=>"mv blib/lib/* $prefix/legacy/");

# Digest::MD5
`mkdir -p PhyloSift/legacy/md5lib`;
add_package(url=>"http://search.cpan.org/CPAN/authors/id/G/GA/GAAS/Digest-MD5-2.52.tar.gz", mv_cmd=>"mv blib/lib/* $prefix/legacy/md5lib/");


# support for App::Cmd
add_package(url=>"http://search.cpan.org/CPAN/authors/id/A/AD/ADAMK/List-MoreUtils-0.33.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DO/DOY/Package-Stash-0.35.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/Data-OptList-0.108.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/Getopt-Long-Descriptive-0.092.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/String-RewritePrefix-0.006.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/App-Cmd-0.318.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DR/DROLSKY/Class-Load-0.20.tar.gz");

add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/Sub-Install-0.926.tar.gz", mv_cmd=>"mv blib/lib/Sub $prefix/lib");

add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DR/DROLSKY/Module-Implementation-0.06.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/Z/ZE/ZEFRAM/Module-Runtime-0.013.tar.gz", mv_cmd=>"mv blib/lib/Module/* $prefix/lib/Module/");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DO/DOY/Try-Tiny-0.16.tar.gz");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/S/SI/SIMONW/Module-Pluggable-4.3.tar.gz", mv_cmd=>"mv blib/lib/Module/* $prefix/lib/Module/; mv blib/lib/Devel $prefix/lib/");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DR/DROLSKY/Package-DeprecationManager-0.13.tar.gz", mv_cmd=>"mv blib/lib/Package/* $prefix/lib/Package");
add_package(url=>"http://search.cpan.org/CPAN/authors/id/R/RJ/RJBS/Sub-Exporter-0.984.tar.gz", mv_cmd=>"mv blib/lib/Sub/* $prefix/lib/Sub");

add_package(url=>"http://search.cpan.org/CPAN/authors/id/A/AD/ADAMK/Params-Util-1.07.tar.gz");

# requires Build.pl
# provides a pure perl impl
add_package(url=>"http://search.cpan.org/CPAN/authors/id/D/DR/DROLSKY/Params-Validate-1.06.tar.gz", mv_cmd=>"mv blib/lib/Params/Validate* $prefix/lib/Params/", build_pl=>1);

# package everything up and datestamp it
my @timerval = localtime();
my $datestr  = ( 1900 + $timerval[5] );
$datestr .= 0 if $timerval[4] < 9;
$datestr .= ( $timerval[4] + 1 );
$datestr .= 0 if $timerval[3] <= 9;
$datestr .= $timerval[3];
`mv PhyloSift phylosift_$datestr`;
`tar cjf phylosift_$datestr.tar.bz2 phylosift_$datestr`;
`rm -rf phylosift_$datestr`;
`echo "phylosift_$datestr" > psversion`;
exit 0;

sub add_package {
	my %args = @_;
	my $url = $args{url};
	my $make_opts = $args{make_opts} || "";
	my $mv_cmd = $args{mv_cmd};
	my $build_pl = $args{build_pl} || 0;
	my $fname = $url;
	$fname =~ s/http.+\///g;
	my $bname = basename($fname, ".tar.gz");
	
	`curl -LO $url`;
	`tar xzf $fname`;
	chdir($bname);
	if($build_pl){
		`perl Build.PL`;
		`./Build`;
		`$mv_cmd`;
	}else{
		`perl Makefile.PL $make_opts`;
		`make`;
		if(defined $mv_cmd){
			`$mv_cmd`;
		}else{
			`mv -f -u blib/lib/* $prefix/lib/`;
		}
	}
	chdir("..");
}
