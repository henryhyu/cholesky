#!/bin/bash

out_name="out.json"
iter=5
for dir in $@
do
	# fills outfile
	for mtx in ./"$dir"*.mtx
	do
		./solve "$mtx" $dir$out_name $iter
	done
	# creates bargraph .png's
	python bargraph.py $dir $dir$out_name
#	rm $dir$out_name
done
