
for dir in ../TestData/TplusE
do
	for mtx in ./"$dir"/*.mtx
	do
		echo "$mtx"
		./solve "$mtx"
	done
done