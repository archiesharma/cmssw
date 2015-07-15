#!/bin/bash

for s in "5" "10" "50" "100" "500" "1000"; do

	for k in "Case1" "Case2" "Case3" "Case4"; do

		name1="MuMinus"$s"_"$k"_ReRecoPU.txt"
		name2="MuPlus"$s"_"$k"_ReRecoPU.txt"
		name3="MuMinus"$s"_"$k"_ReRecoPU_tot.txt"

		cat $name1 $name2 > ./sum/$name3

	done

done
