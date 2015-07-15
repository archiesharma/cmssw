#!/bin/bash

for i in `ls -l | awk '{print $9}' | grep 'Mu'`; do

	echo $i
	crab -submit -c $i
	#crab -kill all -c $i
	#crab -status -c $i
	#crab -get -c $i
	#crab -publish -c $i

done
