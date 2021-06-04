import numpy as np
import sys

file = sys.argv[1]

d = []
fh = open(file, 'r')
for line in fh:
	line = line.strip()
	d.append(int(line.split('\t')[2]))
result = round(np.mean(d),3)
print ('Mean coverage = ', result)
