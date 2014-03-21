package Phylosift::BeastInterface;
use warnings;
use strict;
use Cwd;
use Bio::SeqIO;
use Bio::AlignIO;

our $VERSION = "v1.0.1";

=head1 NAME

Phylosift::BeastInterface - functions to interface with the BEAST de novo metagenome phylogeny module

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Phylosift::Phylosift;

    my $foo = Phylosift::Phylosift->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 export
    
    Args : $outputFile 

    Writes BEAST model XML to the specified file
=cut

sub Export($$$) {
	my %args       = @_;
	my $self       = $args{self} || miss("self");
	my $markRef    = $args{marker_reference} || miss("marker_reference");
	my $outputFile = $args{output_file} || miss("output_file");
	foreach my $marker ( @{$markRef} ) {
		my $trimfinalFastaFile = Phylosift::Utilities::get_trimfinal_fasta_marker_file( self => $self, marker => $marker );
		my $trimfinalFile = Phylosift::Utilities::get_trimfinal_marker_file( self => $self, marker => $marker );
		my $treeFile = Phylosift::Utilities::get_tree_marker_file( self => $self, marker => $marker );
		my $treeStatsFile = Phylosift::Utilities::get_tree_stats_marker_file( self => $self, marker => $marker );
		my $readAlignmentFile = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker );
		if ( !-e "$trimfinalFastaFile" ) {
			`cp "$trimfinalFile" "$trimfinalFastaFile"`;
		}
	}

	# read all alignment input files
	# record which taxa are in each file
	my @alignio;
	my @metaalignio;
	my $inputcount = scalar(@ARGV);
	my %refseqs;
	my @metareadseqs;
	foreach my $marker ( @{$markRef} ) {
		my $trimfinalFastaFile = Phylosift::Utilities::get_trimfinal_fasta_marker_file( self => $self, marker => $marker );
		my $readAlignmentFile = $self->{"alignDir"}."/".Phylosift::Utilities::get_aligner_output_fasta_AA( marker => $marker );

		# read the alignment of reads
		my $in = Bio::AlignIO->new( -file => $readAlignmentFile, '-format' => 'fasta' );
		while ( my $aln = $in->next_aln() ) {
			$aln->sort_alphabetically();
			push( @metaalignio, $aln );
			foreach my $seq ( $aln->each_seq() ) {
				push( @metareadseqs, $seq->id() );
			}
		}

		# read the reference sequence alignment
		my $inphylosift = Bio::AlignIO->new( -file => $trimfinalFastaFile, '-format' => 'fasta' );
		while ( my $aln = $inphylosift->next_aln() ) {
			$aln->sort_alphabetically();
			push( @alignio, $aln );
			foreach my $seq ( $aln->each_seq() ) {
				my $taxon = $seq->desc();
				if ( $taxon =~ /\[Chromosome/ ) {
					$taxon =~ s/^.+?{//g;
					$taxon =~ s/}.+$//g;
				} else {
					$taxon =~ s/.+?\[(.+)\]$/$1/g;
				}
				$taxon =~ s/-//g;
				$refseqs{$taxon} = 0 unless defined( $refseqs{$taxon} );
				$refseqs{$taxon}++;
				$seq->display_id($taxon);
			}
		}
	}
	return if ( @metareadseqs == 0 );    # nothing to see here...move along
	writeXML(
			  xml_file                  => $outputFile,
			  marker_reference          => $markRef,
			  ref_seq_reference         => \%refseqs,
			  meta_Reads_seqs_reference => \@metareadseqs,
			  alignoi_reference         => \@alignio,
			  meta_alignio_reference    => \@metaalignio,
			  cmult                     => 1
	);
}

=head2 writeXML



=cut

sub writeXML {
	my %args            = @_;
	my $xmlfile         = $args{xml_file} || miss("xml_file");
	my $markRef         = $args{marker_reference} || miss("marker_reference");
	my $refseqsref      = $args{ref_seq_reference} || miss("ref_seq_reference");
	my $metareadseqsref = $args{meta_reads_seqs_reference} || miss("meta_reads_seqs_reference");
	my $alignioref      = $args{alignio_reference} || miss("alignio_reference");
	my $metaalignioref  = $args{meta_alignio_reference} || miss("meta_alignio_reference");
	my $cmult           = $args{cmult} || miss("cmult");                                           # 1 for amino acid, 3 for codons
	my %refseqs         = %$refseqsref;
	my @metareadseqs    = @$metareadseqsref;
	my @alignio         = @$alignioref;
	my @metaalignio     = @$metaalignioref;
	my $ALLXML          = ps_open(">$xmlfile");
	print $ALLXML "<beast>\n";
	print $ALLXML "\t<taxa id=\"taxa\">\n";

	foreach my $refseq ( keys(%refseqs) ) {
		print $ALLXML "\t\t<taxon id=\"$refseq\"/>\n";
	}
	foreach my $metaseq (@metareadseqs) {
		print $ALLXML "\t\t<taxon id=\"$metaseq\"/>\n";
	}
	print $ALLXML "\t</taxa>\n";
	print $ALLXML "\t<taxa id=\"referenceTaxa\">\n";
	foreach my $refseq ( keys(%refseqs) ) {
		print $ALLXML "\t\t<taxon idref=\"$refseq\"/>\n";
	}
	print $ALLXML "\t</taxa>\n";
	my @genenames;
	my @markers = @{$markRef};
	for ( my $i = 0; $i < @markers; $i++ ) {
		my $genename = $markers[$i];
		$genename =~ s/\..+//g;
		push( @genenames, $genename );
		print $ALLXML "\t<alignment id=\"$genename.alignment\" dataType = \"amino acid\">\n" if $cmult == 1;
		print $ALLXML "\t<alignment id=\"$genename.alignment\" dataType = \"nucleotide\">\n" if $cmult == 3;

		# first print ref seq alignment
		my $curlen = 0;
		my %seentaxa;
		foreach my $seq ( $alignio[$i]->each_seq() ) {
			my $taxon = $seq->id();
			$taxon =~ s/.+\-//g;
			print $ALLXML "\t\t<sequence>\n";
			print $ALLXML "\t\t<taxon idref=\"$taxon\"/>\n";
			print $ALLXML "\t\t".$seq->seq()."\n";
			print $ALLXML "\t\t</sequence>\n";
			$seentaxa{$taxon} = 1;
			$curlen = length( $seq->seq() );
		}
		foreach my $taxon ( keys(%refseqs) ) {
			next if defined( $seentaxa{$taxon} );

			#			print STDERR "Missing sequence for taxon $taxon in gene $genename.  Filling with gaps.\n";
			print $ALLXML "\t\t<sequence>\n";
			print $ALLXML "\t\t<taxon idref=\"$taxon\"/>\n";
			print $ALLXML "\t\t"."-" x $curlen."\n";
			print $ALLXML "\t\t</sequence>\n";
		}

		# then read seq alignment
		foreach my $seq ( $metaalignio[$i]->each_seq() ) {
			print $ALLXML "\t\t<sequence>\n";
			my $seqid = $seq->id();
			print $ALLXML "\t\t<taxon idref=\"$seqid\"/>\n";
			print $ALLXML "\t\t".$seq->seq()."\n";
			print $ALLXML "\t\t</sequence>\n";
		}
		print $ALLXML "\t</alignment>\n";
	}
	my $alnsuffix = "alignment";
	if ( $cmult == 3 ) {
		foreach my $gene (@genenames) {
			$alnsuffix = "codonAlignment";
			print $ALLXML qq(
			<convert id="$gene.codonAlignment" dataType="codon">
			<alignment idref="$gene.alignment"/>
			</convert>
		);
		}
		print $ALLXML qq(
	<mergePatterns id="allpatterns">
	);
		foreach my $gene (@genenames) {
			print $ALLXML qq(
			<alignment idref="$gene.codonAlignment"/>);
		}
		print $ALLXML qq(
	</mergePatterns>
	);
	}
	foreach my $gene (@genenames) {
		print $ALLXML qq(
		<metagenomeData id="$gene.metagenomeData">
			<taxa idref="referenceTaxa"/>
			<alignment idref="$gene.$alnsuffix"/>
		</metagenomeData>

		<hiddenLinkageModel id="$gene.hiddenStrains" linkageGroupCount="50">
			<metagenomeData idref="$gene.metagenomeData"/>
		</hiddenLinkageModel>
		);
	}

	# this has to get tailored to use one of the hiddenLinkageModels
	print $ALLXML qq(
		<!-- This is a simple constant population size coalescent model              -->
		<!-- that is used to generate an initial tree for the chain.                 -->
		<constantSize id="initialDemo" units="substitutions">
			<populationSize>
				<parameter id="initialDemo.popSize" value="100.0"/>
			</populationSize>
		</constantSize>

		<!-- Generate a random starting tree under the coalescent process            -->
		<!-- No calibration                                                          -->
		<coalescentTree id="startingTree" rootHeight="0.8">
			<hiddenLinkageModel idref="$genenames[0].hiddenStrains"/>
			<constantSize idref="initialDemo"/>
		</coalescentTree>
		<!-- Generate a tree model                                                   -->
		<treeModel id="treeModel">
			<coalescentTree idref="startingTree"/>
			<rootHeight>
				<parameter id="treeModel.rootHeight"/>
			</rootHeight>
			<nodeHeights internalNodes="true">
				<parameter id="treeModel.internalNodeHeights"/>
			</nodeHeights>
			<nodeHeights internalNodes="true" rootNode="true">
				<parameter id="treeModel.allInternalNodeHeights"/>
			</nodeHeights>
		</treeModel>
	);

	# each gene should have its own branch rates
	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(
		<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->
		<discretizedBranchRates id="$gene.branchRates">
			<treeModel idref="treeModel"/>
			<distribution>
				<logNormalDistributionModel meanInRealSpace="true">
					<mean>
						<parameter id="$gene.ucld.mean" value="1.0" lower="0.0" upper="Infinity"/>
					</mean>
					<stdev>
						<parameter id="$gene.ucld.stdev" value="0.1" lower="0.0" upper="Infinity"/>
					</stdev>
				</logNormalDistributionModel>
			</distribution>
			<rateCategories>
				<parameter id="$gene.branchRates.categories" dimension="3346"/>
			</rateCategories>
		</discretizedBranchRates>
		<rateStatistic id="$gene.meanRate" name="$gene.meanRate" mode="mean" internal="true" external="true">
			<treeModel idref="treeModel"/>
			<discretizedBranchRates idref="$gene.branchRates"/>
		</rateStatistic>
		<rateStatistic id="$gene.coefficientOfVariation" name="$gene.coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">
			<treeModel idref="treeModel"/>
			<discretizedBranchRates idref="$gene.branchRates"/>
		</rateStatistic>
		<rateCovarianceStatistic id="$gene.covariance" name="$gene.covariance">
			<treeModel idref="treeModel"/>
			<discretizedBranchRates idref="$gene.branchRates"/>
		</rateCovarianceStatistic>
		);
	}
	print $ALLXML qq(

		<yuleModel id="yule" units="substitutions">
			<birthRate>
				<parameter id="yule.birthRate" value="1.0" lower="0.0" upper="1000000.0"/>
			</birthRate>
		</yuleModel>

		<speciationLikelihood id="speciation">
			<model>
				<yuleModel idref="yule"/>
			</model>
			<speciesTree>
				<treeModel idref="treeModel"/>
			</speciesTree>
		</speciationLikelihood>
	);
	if ( $cmult == 1 ) {
		print $ALLXML qq(
			<!-- The substitution model                                         -->
			<aminoAcidModel id="aa" type="JTT"/>
			<!-- site model                                                              -->
			<siteModel id="siteModel">
				<substitutionModel>
					<aminoAcidModel idref="aa"/>
				</substitutionModel>
			</siteModel>
		);
	} else {
		print $ALLXML qq(
		     <yangCodonModel id="yang" geneticCode="universal">
			<omega>
			    <parameter id="omega" value="1.414"/>
			</omega>
			<kappa>
			    <parameter id="kappa" value="1.234567"/>
			</kappa>
			<frequencyModel id="freqmodel" dataType="codon-universal">
			    <alignment idref="allpatterns"/>
			    <frequencies>
				<parameter id="codonfrequencies" dimension="61"/>
			    </frequencies>
			</frequencyModel>
		    </yangCodonModel>

			<!-- site model                                                              -->
			<siteModel id="siteModel">
				<substitutionModel>
					<yangCodonModel idref="yang"/>
				</substitutionModel>
				<mutationRate>
					<parameter id="siteModel1.mu" value="1.0" lower="0.0"/>
				</mutationRate>
			</siteModel>
		);
	}

	# each gene needs one of these
	foreach my $gene (@genenames) {
		print $ALLXML qq(
		<treeLikelihood id="$gene.treeLikelihood" useAmbiguities="false">
			<hiddenLinkageModel idref="$gene.hiddenStrains"/>
			<treeModel idref="treeModel"/>
			<siteModel idref="siteModel"/>
		);
		print $ALLXML "	<discretizedBranchRates idref=\"allgenes.branchRates\"/>\n";
		print $ALLXML qq(
		</treeLikelihood>
		);
	}
	print $ALLXML qq(	<!-- Define operators                                                        -->
		<operators id="operators">
	);

	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(
			<scaleOperator scaleFactor="0.75" weight="3">
				<parameter idref="$gene.ucld.mean"/>
			</scaleOperator>
			<scaleOperator scaleFactor="0.75" weight="3">
				<parameter idref="$gene.ucld.stdev"/>
			</scaleOperator>
			<upDownOperator scaleFactor="0.75" weight="3">
				<up>
					<parameter idref="$gene.ucld.mean"/>
				</up>
				<down>
					<parameter idref="treeModel.allInternalNodeHeights"/>
				</down>
			</upDownOperator>
			<swapOperator size="1" weight="10" autoOptimize="false">
				<parameter idref="$gene.branchRates.categories"/>
			</swapOperator>
			<randomWalkOperator windowSize="1" weight="10">
				<parameter idref="$gene.branchRates.categories"/>
			</randomWalkOperator>
			<uniformIntegerOperator weight="10">
				<parameter idref="$gene.branchRates.categories"/>
			</uniformIntegerOperator>
		);
	}
	print $ALLXML qq(		<subtreeSlide size="0.08" gaussian="true" weight="15">
				<treeModel idref="treeModel"/>
			</subtreeSlide>
			<narrowExchange weight="15">
				<treeModel idref="treeModel"/>
			</narrowExchange>
			<wideExchange weight="3">
				<treeModel idref="treeModel"/>
			</wideExchange>
			<wilsonBalding weight="3">
				<treeModel idref="treeModel"/>
			</wilsonBalding>
			<scaleOperator scaleFactor="0.75" weight="3">
				<parameter idref="treeModel.rootHeight"/>
			</scaleOperator>
			<uniformOperator weight="30">
				<parameter idref="treeModel.internalNodeHeights"/>
			</uniformOperator>
			<scaleOperator scaleFactor="0.75" weight="3">
				<parameter idref="yule.birthRate"/>
			</scaleOperator>
	);
	foreach my $gene (@genenames) {
		print $ALLXML qq(
			<moveLinkageGroup weight="15.0">
			    <hiddenLinkageModel idref="$gene.hiddenStrains"/>
			</moveLinkageGroup>
			<linkageGroupSwap weight="15.0">
			    <hiddenLinkageModel idref="$gene.hiddenStrains"/>
			</linkageGroupSwap>
		);
	}
	if ( $cmult == 3 ) {
		print $ALLXML qq(
			<scaleOperator scaleFactor="0.75" weight="1.0">
				<parameter idref="siteModel1.mu"/>
			</scaleOperator>
		);
	}
	print $ALLXML qq(	</operators>
		<!-- Define MCMC                                                             -->
		<mcmc id="mcmc" chainLength="10000000" autoOptimize="true">
			<posterior id="posterior">
				<prior id="prior">
	);

	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(
					<gammaPrior shape="0.0010" scale="1000.0" offset="0.0">
						<parameter idref="$gene.ucld.mean"/>
					</gammaPrior>
		);
	}
	print $ALLXML qq(				<speciationLikelihood idref="speciation"/>
				</prior>
				<compoundLikelihood id="likelihood" threads="8">
	);
	foreach my $gene (@genenames) {
		print $ALLXML qq(
					<treeLikelihood idref="$gene.treeLikelihood"/>
		);
	}
	print $ALLXML qq(			</compoundLikelihood>
			</posterior>
			<operators idref="operators"/>
			<!-- write log to screen                                                     -->
			<log id="screenLog" logEvery="1000">
				<column label="Posterior" dp="4" width="12">
					<posterior idref="posterior"/>
				</column>
				<column label="Prior" dp="4" width="12">
					<prior idref="prior"/>
				</column>
				<column label="Likelihood" dp="4" width="12">
					<likelihood idref="likelihood"/>
				</column>
				<column label="rootHeight" sf="6" width="12">
					<parameter idref="treeModel.rootHeight"/>
				</column>
	);

	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(			<column label="$gene.ucld.mean" sf="6" width="12">
					<parameter idref="$gene.ucld.mean"/>
				</column>
		);
	}
	print $ALLXML qq(		</log>
			<!-- write log to file                                                       -->
			<log id="fileLog" logEvery="1000" fileName="frr.hits.trim50.length100.dups.resolved.named.log">
				<posterior idref="posterior"/>
				<prior idref="prior"/>
				<likelihood idref="likelihood"/>
				<parameter idref="treeModel.rootHeight"/>
				<parameter idref="yule.birthRate"/>
	);

	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(			<parameter idref="$gene.ucld.mean"/>
				<parameter idref="$gene.ucld.stdev"/>
				<rateStatistic idref="$gene.meanRate"/>
				<rateStatistic idref="$gene.coefficientOfVariation"/>
				<rateCovarianceStatistic idref="$gene.covariance"/>
		);
	}
	foreach my $gene (@genenames) {
		print $ALLXML qq(			<treeLikelihood idref="$gene.treeLikelihood"/>
		);
	}
	print $ALLXML qq(			<speciationLikelihood idref="speciation"/>
			</log>
			<!-- write tree log to file                                                  -->
			<logTree id="treeFileLog" logEvery="1000" nexusFormat="true" fileName="frr.hits.trim50.length100.dups.resolved.named.trees" sortTranslationTable="true">
				<treeModel idref="treeModel"/>
	);

	#foreach my $gene(@genenames)
	{
		my $gene = "allgenes";
		print $ALLXML qq(			<discretizedBranchRates idref="$gene.branchRates"/>
		);
	}
	print $ALLXML qq(			<posterior idref="posterior"/>
			</logTree>
	);
	foreach my $gene (@genenames) {
		print $ALLXML qq(
			<logHiddenLinkage id="$gene.hiddenLinkageFileLog" logEvery="1000" fileName="linkageGroups.$gene">
				<hiddenLinkageModel idref="$gene.hiddenStrains"/>
			</logHiddenLinkage>
			<logHiddenLinkageTree id="$gene.hiddenLinkageTreeFileLog" logEvery="1000" nexusFormat="true" fileName="hiddenTrees.$gene"
				 sortTranslationTable="true">
			    <treeModel idref="treeModel"/>
			    <posterior idref="posterior"/>
			    <hiddenLinkageModel idref="$gene.hiddenStrains"/>
			</logHiddenLinkageTree>
		);
	}
	print $ALLXML qq(
		</mcmc>
		<report>
			<property name="timer">
				<mcmc idref="mcmc"/>
			</property>
		</report>
	</beast>

	);
	close $ALLXML;
}
