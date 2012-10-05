#! /bin/bash
# Author: Eric Lowe
prefix='/home/elowe/sim_data1' # variable to set path prefix for data files
kounif='knockoutunif_ill_fastq-reads' # variable to set file name for knockoutunif files
koexp='knockoutexp_ill_fastq-reads' # variable to set file name for knockoutexp files
out="$prefix/benchmark_out.txt" # output file
mg='mg' # string comparison variable, for MEGAN
mc='mc' # string comparison variable, for MetaCV
pb='pb' # string comparison variable, for PhymmBL
# declare associative array path
declare -A path
path=([psunif]="$prefix/$kounif.fastq" [mgunif]="$prefix/${kounif}_megan.fastq" [psexp]="$prefix/$koexp.fastq" [mgexp]="$prefix/${koexp}_megan.fastq" [mcunif]="$prefix/${kounif}_mcv.fastq" [mcexp]="$prefix/${koexp}_mcv.fastq" [pbunif]="$prefix/${kounif}_phymmbl.fastq" [pbexp]="$prefix/${koexp}_phymmbl.fastq")
file=(psunif mgunif mcunif psexp mgexp mcexp pbunif pbexp) # array of file "names"
ty=(tophit tophit.recall) # array of file types

if [ -e $out ] ; # does output file exist
then
    echo '' > $out # overwrite contents with blank line
else
    touch $out # else create file
fi

for i in ${file[*]} ; # for each element in file array
do 
    data=${i##[pm][sgcb]} # strips off ps, pb, mc or mg, assigning data as unif or exp
    prog=${i%%[ue]*} # strips off unif or exp, assigning prog as ps, pb, mc or mg
    
    if [ $prog == $mg ] ; # if $prog is mg
    then
	prog='megan' # assign megan to $prog
    elif [ $prog == $mc ] ; # if $prog is mc
    then
	prog='metacv' # assign metacv to $prog
    elif [ $prog == $pb ] ; # if $prog is pb
    then
	prog='phymmbl' # assign phymmbl to $prog
    fi

    for j in ${ty[*]} ; # for each element in ty array
    do 
	if [ -e ${path[$i]}.$j.csv ] ; # if appropriate benchmark output file exists
	then
	    echo "${prog^^*}: $data dataset $j output" >> $out # make prog uppercase, print info line
	    cat ${path[$i]}.$j.csv >> $out # append contents of file to output file
	    echo '' >> $out # add blank line to make easier to read
        fi
    done
done

