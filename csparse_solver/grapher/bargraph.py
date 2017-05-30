import sys
import numpy as np
import matplotlib.pyplot as plt

def autolabel(ax, rects, decimals):
	for rect in rects:
		height = rect.get_height()
		if decimals == 3:
			ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
				'%.3f' % float(height),
				ha='center', va='bottom')
		else:
			ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
				'%.d' % int(height),
				ha='center', va='bottom')


def create_bar(data, xlab, ylab, title, xticks, out):
	n = len(data)
	fig, ax = plt.subplots()
	index = np.arange(n)
	bar_width = 0.35
	opacity = 0.4
	rect = plt.bar(index, data, bar_width, 
		alpha=opacity, 
		color='b')

	autolabel(ax, rect, 3)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	plt.xticks(index, xticks)
	plt.tight_layout()
	# plt.show()
	plt.savefig(out)

def create_bi_bar(data1, data2, xlab, ylab, title, xticks, lab1, lab2, out):
	n = len(data1)
	fig, ax = plt.subplots()
	index = np.arange(n)
	bar_width = 0.35
	opacity = 0.4
	rect1 = plt.bar(index, data1, bar_width, 
		alpha=opacity, 
		color='b',
		label=lab1)

	rect2 = plt.bar(index + bar_width, data2, bar_width, 
		alpha=opacity, 
		color='r',
		label=lab2)

	autolabel(ax, rect1, 1)
	autolabel(ax, rect2, 1)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	plt.xticks(index + bar_width / 2, xticks)
	#plt.tight_layout()
	ax.legend((rect1[0], rect2[0]), ('nnz(A)', 'nnz(R)'))
	# plt.show()
	plt.savefig(out)

# n = number of data entries
def out_parse(file, n):
	data = []
	for i in xrange(n):
		data.append([])
	with open(file, "r") as f:
		counter = 0
		headers = f.readline().split(",")
		for line in f:
			line = line.split(",")
			data[0].append(int(line[0][line[0].index("_") + 10:-4]))
			for i in xrange(1, n):
				data[i].append(float(line[i]))

	return headers, data

if __name__ == "__main__":
	# hard coded
	n = 6
	current_dir = sys.argv[1]
	out_file = sys.argv[2]
	headers, data = out_parse(out_file, n + 2)
	for i in xrange(1, n):
		create_bar(data[i], "Num Edges", "Time (s)", headers[i], data[0], current_dir + headers[i]+".png")
	create_bi_bar(data[6], data[7], "Num Edges", "Nonzeroes", "nnz", data[0], "nnz(A)", "nnz(R)", current_dir + "nnz.png")