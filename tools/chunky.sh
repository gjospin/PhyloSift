#! /bin/bash
# Author: Eric Lowe
# This is a bash script that is meant to make a single sequence_taxa.txt file
# from the multiple files produced by chunking with Phylosift

[ $# -ne 1 ] && { echo "usage: chunky.sh [reads file]" ; exit 1;}

orig="PS_temp/$1/sequence_original.txt"
chunks="PS_temp/$1/sequence_taxa"
target="PS_temp/$1/sequence_taxa.txt"
if [ ! -e $target ] ;
then
    touch $target
fi

if [ -s $target ] ;
then
    mv $target $orig
    touch $target
fi

num=`ls PS_temp/$1/sequence_taxa.* | wc -l`
i=2

cat $chunks.1.txt > $target
while [ $i -lt $num ] ;
do
    cat $chunks.$i.txt >> $target
    let i+=1
done

if [ -s $orig ] ;
then
    diff $orig $target
    if [ $? == "1" ] ;
    then
	echo "Woah, something went wrong! Restoring files"
	mv $orig $target
	exit 1
    fi
    rm $orig
fi