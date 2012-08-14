#! /bin/bash
# Author: Eric Lowe
# Usage: testbench.sh [blastfile] [fastafile]
# For now script requires either a .fa file for fasta file or a .fastq
# with a copy converted to .fa in the same directory.  This is because MEGAN
# only takes .fa, where Phylosift takes .fastq and .fa.  This script WILL NOT
# perform the conversion yet, but that may be implemented later for ease of use. 
# Script must be run in an environment with X11 enabled or X11 forwarding enabled 
# for now.  Beware of running with screen.

# Tests for correct number of arguments and prints usage statement before exiting
# if incorrect number of args supplied
[ $# -ne 2 ] && { echo "usage: testbench.sh [blastfile] [fastafile]" ; exit 1; }

# Declaration of file variables needed for I/O and usage
output=`date +%m%d%y`.csv # MEGAN output .csv file with read ids and tax ids
dirname=`basename $2` # PS_temp directory name to look in for megan_to_ps.pl
meganfile=`date +%m%d%y`_out.rma # MEGAN run output script named with date
ext=`echo $2 | sed 's/.*\.//'` # Gets file extension, for checking if .fastq or .fa
test='fastq' # variable for testing file extension
# fasta file to use; note that sim_data1 is where my file is located, can be changed as needed
# .fa file is required in the directory for script to run properly
file=sim_data1/`basename $2 .fastq`.fa  

if [ -e $meganfile ] ; # if MEGAN has been ran today, do less work
then
    echo -e "open file='$meganfile'\ncollapse rank=Species\nselect rank=Species\nexport what=CSV format=readname_taxonid separator=comma file='/home/elowe/$output'\nquit" > /home/elowe/meg_input.txt
else # MEGAN has not been run today, do long run
    echo -e "load gi2taxfile='/home/elowe/Testing/gi_taxid_nucl.bin'\nimport blastfile='$1' fastafile='$file' meganfile='$meganfile'\ncollapse rank=Species\nselect rank=Species\nexport what=CSV format=readname_taxonid separator=comma file='/home/elowe/$output'\nquit" > /home/elowe/meg_input.txt 
fi # end of if statement

if [ $ext == $test ] ; # if file extension is .fastq
then
    # gi_to_taxid lookup file must be supplied for proper usage.  Mine is located in /home/elowe/Testing
    # but file can be anywhere.  Change as needed so that script finds gi_taxid_nucl.bin
    
    # Runs MEGAN from command line
    MEGAN +g < meg_input.txt 
else
    MEGAN +g < meg_input.txt
fi # end of if statement

# Calls megan_to_ps.pl, which is a perl script (included in Phylosift/tools/) that converts
# MEGAN output to be similar to Phylosift output so that Phylosift Benchmark can be used
# to compare MEGAN to PS.  Output is saved to sequence_taxa.txt in appropriate dir in 
# PS_temp. 
megan_to_ps.pl $output > PS_temp/$dirname/sequence_taxa.txt
# removes unnecessary beginning newline  
sed 's/^\n//' PS_temp/$dirname/sequence_taxa.txt > PS_temp/$dirname/sequence_taxa.txt
# Runs Phylosift benchmark
PhyloSift/bin/phylosift benchmark $2 ~/megan_test 

# End of script