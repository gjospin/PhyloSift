#! /bin/bash
# Author: Eric Lowe
out="/home/elowe/sim_data1/benchmark_out.txt"
mg='mg'
# declare associative array path
declare -A path
path=([psunif]="/home/elowe/sim_data1/knockoutunif_ill_fastq-reads.fastq" [mgunif]="/home/elowe/sim_data1/knockoutunif_ill_fastq-reads_megan.fastq" [psexp]="/home/elowe/sim_data1/knockoutexp_ill_fastq-reads.fastq" [mgexp]="/home/elowe/sim_data1/knockoutexp_ill_fastq-reads_megan.fastq")
file=(psunif mgunif psexp mgexp)
ty=(mass tophit tophit.recall)

if [ -e $out ] ; # does output file exist
then
    echo '' > $out # overwrite contents with blank line
else
    touch $out # else create file
fi

for i in ${file[*]} ; do # for each element in file array
    data=${i##[pm][sg]} # strips off ps or mg, assigning data as unif or exp
    prog=${i%%[ue]*} # strips off unif or exp, assigning prog as ps or mg
    
    if [ $prog == $mg ] ; # if $prog is mg
    then
	prog='megan' # assign megan to $prog
    fi
    for j in ${ty[*]} ; do # for each element in ty array
	echo "${prog^^*}: $data dataset $j output" >> $out # make prog uppercase, print info line
	cat ${path[$i]}.$j.csv >> $out # append contents of file to output file
	echo '' >> $out # add blank line to make easier to read
    done
done

