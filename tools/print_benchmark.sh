#! /bin/bash
# Author: Eric Lowe
out="/home/elowe/sim_data1/benchmark_out.txt"
mg='mg' # string comparison variable, for MEGAN
mc='mc' # string comparison variable, for MetaCV
# declare associative array path
declare -A path
path=([psunif]="/home/elowe/sim_data1/knockoutunif_ill_fastq-reads.fastq" [mgunif]="/home/elowe/sim_data1/knockoutunif_ill_fastq-reads_megan.fastq" [psexp]="/home/elowe/sim_data1/knockoutexp_ill_fastq-reads.fastq" [mgexp]="/home/elowe/sim_data1/knockoutexp_ill_fastq-reads_megan.fastq" [mcunif]="/home/elowe/sim_data1/knockoutunif_ill_fastq-reads_mcv.fastq" [mcexp]="/home/elowe/sim_data1/knockoutexp_ill_fastq-reads_mcv.fastq")
file=(psunif mgunif mcunif psexp mgexp mcexp) # array of file "names"
ty=(tophit tophit.recall) # array of file types

if [ -e $out ] ; # does output file exist
then
    echo '' > $out # overwrite contents with blank line
else
    touch $out # else create file
fi

for i in ${file[*]} ; # for each element in file array
do 
    data=${i##[pm][sgc]} # strips off ps, mc or mg, assigning data as unif or exp
    prog=${i%%[ue]*} # strips off unif or exp, assigning prog as ps, mc or mg
    
    if [ $prog == $mg ] ; # if $prog is mg
    then
	prog='megan' # assign megan to $prog
    fi
    
    if [ $prog == $mc ] ; # if $prog is mc
    then
	prog='metacv' # assign metacv to $prog
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

