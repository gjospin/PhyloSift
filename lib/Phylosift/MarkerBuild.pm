package Phylosift::MarkerBuild;
use Cwd;
use Carp;
use Phylosift::Utilities;
use File::Basename;

=head1 NAME

Phylosift::MarkerBuild - build a seed marker package from an existing multiple sequence alignment

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Given a multiple sequence alignment, functions in this package can build a marker package.
Steps involved include masking the alignment, selecting representative sequences for database searches, building a phylogenetic tree, HMM, and pplacer package.

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 build_marker

=cut

sub build_marker {
	my %args = @_;
	my $aln_file     = $args{alignment};
	my $package_name = $args{name};
	my $cutoff       = $args{cutoff};
	
	my $target_dir   = getcwd() . "/$package_name";
	`mkdir $target_dir` unless -e $target_dir;
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	
	my $fasta_file = "$target_dir/$core.fasta";
	my $seq_count = Phylosift::Utilities::unalign_sequences($aln_file, $fasta_file);
	
	my $masked_aln = "$target_dir/$core.masked";
	mask(file=>$aln_file,output=>$masked_aln);
	
	my $hmm_file = "$target_dir/$core.hmm";
	generate_hmm( $masked_aln, $hmm_file );
	
	Phylosift::Utilities::fastpsstockholm($masked_aln, "$target_dir/$core.stk");
	my $stk_aln =  "$target_dir/$core.stk";

	#may need to create an unaligned file for the sequences before aligning them
	my $new_alignment_file = hmmalign_to_model( $hmm_file, $fasta_file, $target_dir, $stk_aln, $seq_count);
	my $clean_aln = "$target_dir/$core.clean";
	my %id_map = mask_and_clean_alignment( $new_alignment_file , $clean_aln );
	my ( $fasttree_file, $tree_log_file ) = generate_fasttree( $clean_aln, $target_dir );
	
	#need to generate representatives using PDA
	my $rep_file = get_representatives_from_tree( $fasttree_file, $target_dir, $cutoff );
	
	#need to read the representatives picked by PDA and generate a representative fasta file
	my $rep_fasta = get_fasta_from_pda_representatives( $rep_file, $target_dir, $fasta_file, \%id_map );
	
	#use taxit to create a new reference package required for running PhyloSift
	#needed are : 1 alignment file, 1 representatives fasta file, 1 hmm profile, 1 tree file, 1 log tree file.
	`cd $target_dir;taxit create -c -d "Creating a reference package for PhyloSift for the $core marker" -l $core -f $clean_aln -t $target_dir/$core.tree -s $target_dir/$core.log -Y FastTree -P $core`;
	`rm $target_dir/$core.pda`;
	`rm $target_dir/$core.tree`;
	`rm $target_dir/$core.log`;
	`rm $target_dir/$core.aln`;
	`rm $target_dir/$core.fasta`;
	`rm $clean_aln`;
	`rm $masked_aln`;
	`mv $target_dir/$core/* $target_dir`;
	`rm -rf $target_dir/$core`;

}


=head2 generate_hmm

input: alignment_file, target_directory
generates a HMM profile from an alignment in FASTA format (arg) using hmmbuild.  The hmm is placed in the target_directory

=cut

sub generate_hmm {
	my $file_name  = shift;
	my $hmm_name = shift;
	`$Phylosift::Utilities::hmmbuild --informat afa $hmm_name $file_name`;
}

=head2 hmmalign_to_model

input : hmm_profile,sequence_file,target_dir
Aligns sequences to an HMM model and outputs an alignment
=cut

sub hmmalign_to_model {
	my $hmm_profile   = shift;
	my $sequence_file = shift;
	my $target_dir    = shift;
	my $ref_ali       = shift;
	my $seq_count     = shift;
	my ( $core_name, $path, $ext ) = fileparse( $sequence_file, qr/\.[^.]*$/ );
	open(ALNOUT, ">$target_dir/$core_name.aln");
	open(ALNIN, "$Phylosift::Utilities::hmmalign --mapali $ref_ali --trim --outformat afa $hmm_profile $sequence_file |");
	my $s = 0;
	while(my $line = <ALNIN>){
		if($line=~/^>/){
			$s++;
			last if $s>$seq_count;
		}
		print ALNOUT $line;
	}
	close ALNOUT;
	return "$target_dir/$core_name.aln";
}

=head2 mask_alignment

input : aln_file
Masks the unaligned columns out of an alignment file. Removes ( and ) from the sequence names 
Also removes duplicate IDs
=cut

sub mask_and_clean_alignment {
	my $aln_file   = shift;
	my $output_file = shift;
	my %id_map;	# will store a map of unique IDs to sequence names

	#    open(FILEIN,$aln_file) or carp("Couldn't open $aln_file for reading \n");
	my $in = Bio::SeqIO->new( -file => $aln_file );
	my %s = ();    #hash remembering the IDs already printed
	open( FILEOUT, ">$output_file" ) or carp("Couldn't open $output_file for writing\n");
	my $seq_counter = 0;
	while ( my $seq_object = $in->next_seq() ) {
		my $seq = $seq_object->seq;
		my $id  = $seq_object->id;
		my $unique_id = sprintf("%09d", $seq_counter++);
		$id_map{$id}=$unique_id;
		$id  =~ s/\(\)//g;     #removes ( and ) from the header lines
		$seq =~ s/[a-z]//g;    # lowercase chars didnt align to model
		$seq =~ s/\.//g;       # shouldnt be any dots
		print FILEOUT ">" . $unique_id . "\n" . $seq . "\n";
	}

	#    close(FILEIN);
	close(FILEOUT);
	return %id_map;
}

=head2 generate_fasttree

input: alignment_file,target_directory
generates a tree using fasttree and write the output along with the log/info files to the target directory.

=cut

sub generate_fasttree {
	my $aln_file   = shift;
	my $target_dir = shift;
	my ( $core, $path, $ext ) = fileparse( $aln_file, qr/\.[^.]*$/ );
	my %type = Phylosift::Utilities::get_sequence_input_type($aln_file);
	if ( $type{seqtype} eq "dna" ) {
		`$Phylosift::Utilities::fasttree -nt -gtr -log $target_dir/$core.log $aln_file > $target_dir/$core.tree 2> /dev/null`;
	} else {
		`$Phylosift::Utilities::fasttree -gtr -log $target_dir/$core.log $aln_file > $target_dir/$core.tree 2> /dev/null`;
	}
	return ( "$target_dir/$core.tree", "$target_dir/$core.log" );
}

=head2 get_representatives_from_tree

input : tree file in newick format, directory to write the output to, pruning threshold
uses the PDA program to prune a tree to get representative sequences

=cut

sub get_representatives_from_tree {
	my $tree_file  = shift;
	my $target_dir = shift;
	my $cutoff     = shift;
	my ( $core, $path, $ext ) = fileparse( $tree_file, qr/\.[^.]*$/ );

	#get the number of taxa in the tree
	my $taxa_count = 0;
	my $input_tree = new Bio::TreeIO( -file => $tree_file, -format => "newick" );
	while ( my $tree = $input_tree->next_tree ) {
		for my $node ( $tree->get_nodes ) {
			if ( $node->is_Leaf ) {
				$taxa_count++;
			}
		}
	}

	#pda doesn't seem to want to run if $taxa_count is the number of leaves. Decrementing to let pda do the search.
	$taxa_count--;
	my $pda_cmd = "cd $target_dir;$Phylosift::Utilities::pda -g -k $taxa_count -minlen $cutoff $tree_file $target_dir/$core.pda";
	`$pda_cmd`;
	return "$target_dir/$core.pda";
}

=head2 get_fasta_from_pda_representatives 

input pda file and reference fasta file
reads the selected representatives from the pda file and prints the sequences to a new fasta file

=cut

sub get_fasta_from_pda_representatives {
	my $pda_file        = shift;
	my $target_dir      = shift;
	my $reference_fasta = shift;
	my $id_map     = shift;
	my ( $core, $path, $ext ) = fileparse( $pda_file, qr/\.[^.]*$/ );

	#reading the pda file to get the representative IDs
	open( REPSIN, $pda_file ) or carp("Could not open $pda_file\n");
	my $taxa_number   = 0;
	my %selected_taxa = ();
	while (<REPSIN>) {
		chomp($_);
		if ( $_ =~ m/optimal PD set has (\d+) taxa:/ ) {
			$taxa_number = $1;
		} elsif ( $_ =~ m/Corresponding sub-tree:/ ) {
			last;
		} elsif ( $taxa_number != 0 && scalar( keys(%selected_taxa) ) < $taxa_number ) {
			$_ =~ m/^(\S+)$/;
			$selected_taxa{$1} = 1;
		}
	}
	close(REPSIN);	

	#reading the reference sequences and printing the selected representatives using BioPerl
	my $reference_seqs        = Bio::SeqIO->new( -file => $reference_fasta,         -format => "FASTA" );
	my $representatives_fasta = Bio::SeqIO->new( -file => ">$target_dir/$core.rep", -format => "FASTA" );
	while ( my $ref_seq = $reference_seqs->next_seq ) {
		if ( exists $selected_taxa{ $id_map->{$ref_seq->id} } ) {
			$representatives_fasta->write_seq($ref_seq);
		}
	}
	return "$target_dir/$core.rep";
}


=head2 mask

Remove columns from a sequence alignment that contain too many gaps or that are otherwise unreliable

=cut

sub mask {
	my %args = @_;
	my $infile = $args{file};
	my $outfile = $args{output};
	
	my $gap_cutoff = 10;
	my $cutoff = 10;
	
	my $maskcont = martin_mask( $args{file}, $cutoff, $gap_cutoff );

	my %maskseq;

	my @ori_order;
	$maskcont =~ s/^>//;
	my @tempmask = split( />/, $maskcont );
	foreach my $tempmask (@tempmask) {
		my @templine = split( /\n/, $tempmask );
		my $this_ID  = shift(@templine);
		my $this_seq = join( '', @templine );
		$this_seq =~ s/\s+//g;
		my @this_seq = split( //, $this_seq );
		if ( $this_ID ne '_mask' ) { push( @ori_order, $this_ID ) }
		$maskseq{$this_ID} = \@this_seq;
	}
	my %trimseq;

	while ( @{ $maskseq{_mask} } ) {
		my $switch = shift( @{ $maskseq{_mask} } );

		foreach my $key ( keys %maskseq ) {
			if ( $key ne '_mask' ) {
				my $aa = shift( @{ $maskseq{$key} } );
				if ( $switch == 1 ) { $trimseq{$key} .= $aa; }
			}
		}
	}

	# writing out a fasta
	open( TRIMOUT, "> $outfile" ) || die "cannot output aligment after trimming\n";
	foreach my $key (@ori_order) {
		print TRIMOUT ">" . $key . "\n";
		my $i;
		for ( $i = 0 ; $i <= length( $trimseq{$key} ) ; $i += 80 ) {
			my $substr = substr( $trimseq{$key}, $i, 80 );
			print TRIMOUT $substr . "\n";
		}
	}
	close TRIMOUT;

}

sub martin_mask {

	my ( $input_file, $cutoff, $opt_g ) = @_;

	if ( !( length($opt_g) > 1 ) ) { $opt_g = 101; }

	my $self = {};
	$self->{inputfile}  = $input_file;
	$self->{cutoff}     = $cutoff;
	$self->{gap_cutoff} = $opt_g;
	$self               = &ReadMatrix($self);
	$self               = &ReadAlignment($self);
	$self               = &CalculateScore($self);
	$self               = &Mask($self);
	my $return_mask = &Output($self);

	return $return_mask;
}

sub ReadMatrix {
	my $self = shift;

	my @aa;
	my %matrix;
	my ( $matrix_name, $i, $j );
	my %min        = ();
	my $ori_matrix = &ori_matrix();
	my @mx         = split( /\n/, $ori_matrix );

	while (@mx) {
		$_ = shift(@mx);
		if (/amino_acid_order = "(.+)"/) {
			@aa = split //, $1;
		} elsif (/MATRIX (.+)\[/) {
			$matrix_name = $1;
			$j           = 0;
		} elsif (/\d,/) {
			s/(\s|\}|;)//g;
			my @score = split /,/;
			for ( $i = 0 ; $i <= $#score ; $i++ ) {
				$matrix{$matrix_name}{"$aa[$i]$aa[$j]"} = $score[$i] * 100;
				$matrix{$matrix_name}{"$aa[$j]$aa[$i]"} = $score[$i] * 100;
			}
			$j++;
		}
	}

	for my $key ( keys %matrix ) {
		for my $pair ( keys %{ $matrix{$key} } ) {
			if ( !$min{$key} ) {
				$min{$key} = $matrix{$key}{$pair};
			} elsif ( $matrix{$key}{$pair} < $min{$key} ) {
				$min{$key} = $matrix{$key}{$pair};
			}
		}
	}

	foreach my $key ( keys %min ) {
		if ( $min{$key} < 0 ) {
			foreach my $pair ( keys %{ $matrix{$key} } ) {
				$matrix{$key}{$pair} -= $min{$key};
			}
		}
	}

	$self->{matrix} = \%matrix;
	$self->{aa}     = \@aa;
	return $self;

}

sub ReadAlignment {

	my $self = shift;
	my %seq;

	my @order;

	my $id;
	my $input_file = $self->{inputfile};
	open( MASKIN, $input_file ) || die "can't open $input_file\n";
	while (<MASKIN>) {
		chop;
		if (/%([\S]+)/) {
			$id = $1;
			push( @order, $id );
		} elsif (/>([\S]+)/) {
			$id = $1;
			push( @order, $id );
		} else {
			$_ =~ s/\s+//g;
			$_ =~ s/\./-/g;	#AED: treat . and ? as gaps
			$_ =~ s/\?/-/g;
			$seq{$id} .= uc($_);
		}
	}
	close(MASKIN);

	$self->{seq}   = \%seq;
	$self->{order} = \@order;
	return $self;

}

sub CalculateScore {

	my $self = shift;
	my %seq  = %{ $self->{seq} };
	my $seqlength;
	my %seqArray;
	my %matrix = %{ $self->{matrix} };
	my %local_score;
	my %column_score;
	my @aa = @{ $self->{aa} };

	for ( keys %seq ) {
		$seqArray{$_} = [ split //, $seq{$_} ];
		$seqlength = length( $seq{$_} );
	}
	if ( $seqlength <= 20 ) {
		carp "Alignment is too short\n";
	}

	my $numseq = scalar keys %seq;

	for my $i ( 0 .. $seqlength - 1 ) {
		my ( %freq, %profile, %dist, %seqVector ) = ();
		my $dist_mean = 0;
		my $number    = 0;
		my $dist_median;
		my $col;
		my $number_gap;

		foreach my $sequence ( keys %seqArray ) {
			$col .= $seqArray{$sequence}[$i];

			if ( $seqArray{$sequence}[$i] ne "-" ) {
				$freq{ $seqArray{$sequence}[$i] }++;
				$number++;
			} else {
				$number_gap++;
			}
		}

		for my $aa1 (@aa) {
			for my $aps (@aa) {
				$profile{$aa1} += $freq{$aps} * $matrix{"gon250mt"}{"$aa1$aps"} if exists $freq{$aps};
			}

			#			$profile{$aa1} /= $number;
			$profile{$aa1} /= $numseq;
		}

		for my $sequence ( keys %seqArray ) {
			my $c = $seqArray{$sequence}[$i];
			if ( $c ne '-' ) {
				for (@aa) {
					$seqVector{$_} = $matrix{"gon250mt"}{"$c$_"};
				}

				for (@aa) {

					my $diff = $profile{$_} - $seqVector{$_};
					$diff /= 1000;

					$dist{$sequence} += $diff * $diff;
				}

				$dist{$sequence} = sqrt( $dist{$sequence} );
				$dist_mean += $dist{$sequence};

				#			$number ++;
			}
		}

		if ($number) {
			$dist_mean /= $number;
		} else {
			$dist_mean = 0;
		}

		my @sort = sort { $a <=> $b } ( values %dist );
		my $t = $number % 2;
		if ( $t == 0 ) {
			$dist_median = ( $sort[ $number / 2 - 1 ] + $sort[ $number / 2 ] ) / 2;
		} else {
			$dist_median = $sort[ ( $number - 1 ) / 2 ];
		}

		#		$column_score2{$i} = exp (-$dist_mean/4) * 100 *$number/$numseq;

		#		$column_score{$i} = exp (-$dist_median/1.82) * 100 * ($number/$numseq);   # 1.82 make total random alignment score 1.0
		$column_score{$i} = exp( -$dist_median / 3 ) * 100 * ( $number / $numseq );
		my $gap_percent = $number_gap / $numseq * 100;
		my $gap_cutoff  = $self->{gap_cutoff};
		if ( ( $gap_cutoff > 0.00001 ) && ( $gap_percent > $gap_cutoff ) ) {
			$self->{rid_gap}->{$i} = 1;
		}

	}

	$self->{column_score} = \%column_score;
	$self->{local_score}  = \%local_score;
##$self->{column_score2}=\%column_score2;
	$self->{seqlength} = $seqlength;
	$self->{seqArray}  = \%seqArray;
	return $self;

}

sub Mask {
	my $self = shift;

	my %matrix       = %{ $self->{matrix} };
	my %seq          = %{ $self->{seq} };
	my %seqArray     = %{ $self->{seqArray} };
	my %column_score = %{ $self->{column_score} };
	my %local_score  = %{ $self->{local_score} };
	my $mask;

	my $seqlength = $self->{seqlength};

	for my $i ( 0 .. $seqlength - 1 ) {

		if ( $i <= 2 or $i >= $seqlength - 3 ) {
			$local_score{$i} = $column_score{$i};
		} else {
			$local_score{$i} =
			  ( $column_score{ $i - 3 } +
				2 * $column_score{ $i - 2 } +
				3 * $column_score{ $i - 1 } +
				4 * $column_score{$i} +
				3 * $column_score{ $i + 1 } +
				2 * $column_score{ $i + 2 } +
				1 * $column_score{ $i + 3 } ) / 16;
		}

		if ( $column_score{$i} == 0 ) { $local_score{$i} = 0; }
		elsif ( $local_score{$i} / $column_score{$i} > 3 ) {
			my $score_l = $column_score{ $i - 3 } + $column_score{ $i - 2 } + $column_score{ $i - 1 };
			my $score_r = $column_score{ $i + 1 } + $column_score{ $i + 2 } + $column_score{ $i + 3 };
			if ( ($score_r) && ( $score_l / $score_r > 20 or $score_l / $score_r < 0.05 ) ) {
				$local_score{$i} = $column_score{$i};
			}
		}
###########################################cutoff #######################################

		if ( ( $local_score{$i} >= $cutoff ) && ( !( $self->{rid_gap}->{$i} ) ) ) {
			$mask .= "1";
		} else {
			$mask .= "0";
		}

		#		print $i+1,"\t";
		#		printf "%0.1f\t%0.1f\n", $column_score{$i},$local_score{$i};
	}

	$self->{mask} = $mask;

	return $self;

}

sub Output {
	my $self = shift;

	my %seq = %{ $self->{seq} };

	my $mask = $self->{mask};

	my @order = @{ $self->{order} };

	my $return_mask = '';

	$seq{"_mask"} = $mask;
	push( @order, '_mask' );

	foreach my $key (@order) {
		$return_mask .= ">$key\n";
		for ( my $i = 0 ; $i < length( $seq{$key} ) ; $i += 60 ) {
			$return_mask .= substr( $seq{$key}, $i, 60 ) . "\n";
		}
	}

	return $return_mask;
}

sub ori_matrix {

	my $matrix = "amino_acid_order = \"ABCDEFGHIKLMNPQRSTVWXYZ\";

MATRIX gon250mt[]={
66,
35,50,
38,35,127,
33,50,13,82,
35,50,15,68,74,
19,35,29,5,9,97,
54,35,21,35,29,0,94,
29,35,26,37,37,34,25,90,
29,35,27,9,17,41,5,20,77,
32,35,16,38,43,13,27,54,21,72,
27,35,25,8,16,48,5,22,69,21,77,
30,35,29,15,21,45,11,26,67,25,69,79,
33,50,23,65,56,14,37,43,16,40,15,20,76,
37,35,14,30,31,9,24,27,17,31,19,19,29,101,
33,50,19,56,62,17,28,43,22,45,24,28,55,33,68,
31,35,20,33,37,13,28,54,19,68,20,23,37,29,45,82,
58,35,35,38,36,16,53,33,23,35,21,25,41,37,36,33,65,
54,35,31,35,34,20,43,33,31,35,26,31,38,35,35,33,60,67,
35,35,35,15,22,35,13,21,71,23,62,61,20,23,25,21,28,35,73,
11,35,28,0,6,74,8,29,23,11,30,28,11,1,17,24,13,11,17,145,
35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,50,
20,35,31,16,17,84,8,49,30,21,35,33,25,14,23,23,22,22,27,78,35,102,
35,50,35,50,50,35,35,35,35,35,35,35,50,35,50,35,35,35,35,35,35,35,50,";
	return $matrix;

}

=head1 AUTHOR

Martin Wu, C<< <> >>
Dongying Wu, C<< <> >>
Aaron Darling, C<< <aarondarling at ucdavis.edu> >>
Guillaume Jospin, C<< <gjospin at ucdavis.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-phylosift-phylosift at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Phylosift-Phylosift>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Phylosift::MarkerBuild


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Phylosift-Phylosift>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Phylosift-Phylosift>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Phylosift-Phylosift>

=item * Search CPAN

L<http://search.cpan.org/dist/Phylosift-Phylosift/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2011 Aaron Darling and Guillaume Jospin.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation.

See http://dev.perl.org/licenses/ for more information.


=cut

1;    # End of Phylosift::MarkerBuild.pm
