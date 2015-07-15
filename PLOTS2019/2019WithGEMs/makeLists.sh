#!/bin/bash

path=$1
root=$(echo `pwd`)
toFind=$2

for s in `ls -l $root/$path | awk '{print $9}' | grep $toFind`; do

	echo $s
	rm $path/$s'.txt'

	for k in `ls -l $path/$s | awk '{print $9}' | grep 'log'`; do

		for i in `cat $path/$s/$k/crab.log | grep 'validationEDM' | grep 'lfn'`; do

			echo $i | grep 'store' >> ./Convert/txt/$s'.txt'

		done

	done

done
