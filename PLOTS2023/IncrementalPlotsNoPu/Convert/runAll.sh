#!/bin/bash

root=$(echo `pwd`)

for i in `ls -l $root/"txt" | grep 'txt' | awk '{print $9}'`; do

	echo $i
	i2=`echo $i | sed "s/.txt$//"`
	echo "_"$i2
	cmsRun EDMtoMEConverter_cfg.py fileList=$root/"txt"/$i output=`echo "_"$i2`

done
