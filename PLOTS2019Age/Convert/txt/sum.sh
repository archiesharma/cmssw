#!/bin/bash

for s in "5" "10" "50" "100" "500" "1000"; do

	for k in "Official"; do

		name1="MuMinus"$s"_"$k"_2019Age.txt"
		name2="MuPlus"$s"_"$k"_2019Age.txt"
		name3="MuMinus"$s"_"$k"_2019Age_tot.txt"

		cat $name1 $name2 > ./sum/$name3

	done

done
