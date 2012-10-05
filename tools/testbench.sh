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

# Tests for X Windowing system since MEGAN requires this to run properly
[ "$DISPLAY" ] || { echo "Need to enable X Windowing!" ; exit 1; }

# Declaration of file variables needed for I/O and usage
dirname=`basename $2` # PS_temp directory name to look in for megan_to_ps.pl
output=$dirname.`date +%m%d%y`.csv # MEGAN output .csv file with read ids and tax ids
meganfile=$dirname.`date +%m%d%y`_out.rma # MEGAN run output script named with date
ext=`echo $2 | sed 's/.*\.//'` # Gets file extension, for checking if .fastq or .fa
test='fastq' # variable for testing file extension
# fasta file to use; note that sim_data1 is where my file is located, can be changed as needed
# .fa file is required in the directory for script to run properly

if [ -e $meganfile ] ; # if MEGAN has been ran today, do less work
then # rank is given on command line, first letter MUST be capital; i.e. Species or Subspecies
echo -e "open file='$meganfile'\ncollapse rank=Species\nselect nodes=all\nexport what=CSV format=readname_taxonid separator=comma file='/home/elowe/$output'\nquit" > /home/elowe/meg_input.txt
else # MEGAN has not been run today, do long run
    if [ $ext == $test ] ; # if file extension is .fastq
    then
        file=sim_data1/`basename $2 .fastq`.fa
        # gi_to_taxid lookup file must be supplied for proper usage.  Mine is located in /home/elowe/Testing
        # but file can be anywhere.  Change as needed so that script finds gi_taxid_nucl.bin
        echo -e "load gi2taxfile='/home/elowe/Testing/gi_taxid_nucl.bin'\nimport blastfile='$1' fastafile='$file' meganfile='$meganfile'\ncollapse rank=Species\nselect nodes=all\nexport what=CSV format=readname_taxonid separator=comma file='/home/elowe/$output'\nquit" > /home/elowe/meg_input.txt
    else
        file=$dirname
        echo -e "load gi2taxfile='/home/elowe/Testing/gi_taxid_nucl.bin'\nimport blastfile='$1' fastafile='$file' meganfile='$meganfile'\ncollapse rank=Species\nselect nodes=all\nexport what=CSV format=readname_taxonid separator=comma file='/home/elowe/$output'\nquit" > /home/elowe/meg_input.txt
    fi # end of if file extension is .fastq   
fi # end of if MEGAN has been run today

# Runs MEGAN from command line
MEGAN +g < meg_input.txt
# Calls megan_to_ps.pl, which is a perl script (included in Phylosift/tools/) that converts
# MEGAN output to be similar to Phylosift output so that Phylosift Benchmark can be used
# to compare MEGAN to PS.  Output is saved to sequence_taxa.txt in appropriate dir in 
# PS_temp. 
if [ -e PS_temp/$dirname/sequence_taxa.txt ] ; # if file exists
then
    megan_to_ps.pl $output > PS_temp/$dirname/sequence_taxa.txt
else # make file so that script can run
    touch PS_temp/$dirname/sequence_taxa.txt
    megan_to_ps.pl $output > PS_temp/$dirname/sequence_taxa.txt
fi

# removes unnecessary beginning newline  
sed -i 's/^\n//' PS_temp/$dirname/sequence_taxa.txt
# Runs Phylosift benchmark
if [ -s PS_temp/$dirname/sequence_taxa.txt ] ; # if file exists and is not empty
then
    phylosift benchmark $2 ~/megan_test 
else # file exists but is empty
    echo "Something went wrong!"
    exit 1
fi
# End of script
