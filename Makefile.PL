use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Amphora2::Amphora2',
    AUTHOR              => q{Aaron Darling <aarondarling@ucdavis.edu> and Guillaume Jospin <gjospin@ucdavis.edu>},
    VERSION_FROM        => 'lib/Amphora2/Amphora2.pm',
    ABSTRACT_FROM       => 'lib/Amphora2/Amphora2.pm',
    LICENSE		=> 'GPL',
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Amphora2-Amphora2-*' },
);
