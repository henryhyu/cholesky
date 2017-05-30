#!/bin/bash

out_name="out.csv"
for dir in $@
do
	# header
	echo "filename,Ordering,NumericFactorization,TriangularSolve,TotalTime,nnzRatio,nnz(A),nnz(R)" > $dir$out_name
	# fills outfile
	for mtx in ./"$dir"*.mtx
	do
		./solve "$mtx" $dir$out_name
	done
	# creates bargraph .png's
	python bargraph.py $dir $dir$out_name
#	rm $dir$out_name
done
