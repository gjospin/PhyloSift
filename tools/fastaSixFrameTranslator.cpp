/**
 * Program to read a FastA file
 * and perform 6-frame translation of paths embedded in the graph
 * @author Aaron Darling
 */

#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnFastTranslator.h"

using namespace std;
using namespace genome;


int main(int argc, char* argv[])
{
	if(argc < 2){
		cerr << "Usage: fastaSixFrameTranslator <FastA file>\n";
		return -1;
	}
	gnSequence gns;
	try{
		gns.LoadSource( argv[1] );
	}catch(...){
		cerr << "Unable to read FastA file \"" << argv[1] << "\"\n";
		return -1;
	}
	for( int cI = 0; cI < gns.contigListSize(); cI++ ){
		string dna = gns.contig(cI).ToString();
		string revdna = dna;
		gnFilter::DNAComplementFilter()->ReverseFilter(revdna);
		for( int i=0; i<3; i++ ){
			string fwd = dna.substr(i);
			string rev = revdna.substr(i);
			gnFastTranslator::DNAProteinTranslator()->Filter(fwd);
			gnFastTranslator::DNAProteinTranslator()->Filter(rev);
			for(size_t j=0; j<fwd.size(); j++)
				if(fwd[j]=='.')	fwd[j] = 'X';
			for(size_t j=0; j<rev.size(); j++)
				if(rev[j]=='.')	rev[j] = 'X';
			string cname = gns.contigName(cI);
			cout << ">" << cname << "_f" << i << endl;
			cout << fwd << endl;
			cout << ">" << cname << "_r" << i << endl;
			cout << rev << endl;
		}
	}
	return 0;
}

