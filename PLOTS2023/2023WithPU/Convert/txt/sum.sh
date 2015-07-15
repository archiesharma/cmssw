#!/bin/bash

for s in "5" "10" "50" "100" "500" "1000"; do

	for k in "Case5"; do

		name1="MuMinus"$s"_"$k"_Official.txt"
		name2="MuPlus"$s"_"$k"_Official.txt"
		name3="MuMinus"$s"_"$k"_Official_tot.txt"

		cat $name1 $name2 > ./sum/$name3

	done

done
