
for dir in TplusE T
do
	for mtx in ./"$dir"/*.mtx
	do
		echo "$mtx"
		./solve "$mtx" >> output
	done
done