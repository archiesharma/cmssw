#!/bin/bash

root=$(echo `pwd`)

for i in `ls -l $root/"txt"/sum | grep '_tot.txt' | awk '{print $9}'`; do

	echo $i
	i2=`echo $i | sed "s/_tot.txt$//"`
	echo "_"$i2
	cmsRun EDMtoMEConverter_cfg.py fileList=$root/"txt"/sum/$i output=`echo "_"$i2`

done
