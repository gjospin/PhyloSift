use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Phylosift::Phylosift',
    AUTHOR              => q{Aaron Darling <aarondarling@ucdavis.edu> and Guillaume Jospin <gjospin@ucdavis.edu>},
    VERSION_FROM        => 'lib/Phylosift/Phylosift.pm',
    ABSTRACT_FROM       => 'lib/Phylosift/Phylosift.pm',
    LICENSE		=> 'GPL',
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Phylosift-Phylosift-*' },
);
