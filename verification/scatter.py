import numpy as numpy
import math
import matplotlib.pyplot as plt
import scipy as sp
import array
import matplotlib.ticker as ticker
from matplotlib import gridspec


plt.rc('text', usetex=True)
font = {'family' : 'sans-serif', 'sans-serif' : 'Helvetica',
          'weight' : 'bold',
          'size'   : 15.}
plt.rc('font', **font )
 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

gs = gridspec.GridSpec(1, 1,
                       #width_ratios=[1,1]
                       )
fig = plt.figure(figsize=(5.6,4.6),dpi=500)

p1 = fig.add_subplot(gs[0])

f1 = open('vector1_initial_clean.txt', 'r')
f2 = open('vector1_results_clean.txt', 'r')

xi=[]
yi=[]

xr=[]
yr=[]

for line in f1:
	xi.append(float(line.split()[0]))
	yi.append(float(line.split()[1]))

for line in f2:
	xr.append(float(line.split()[0]))
	yr.append(float(line.split()[1]))

p1.scatter(xi,yi,marker='o', c='black',label='Initial Eigenvalues',s=70);
p1.scatter(xr,yr, marker='+', c='r',label='Computed Eigenvalues',s=70);

p1.set_xlabel('Real Axis', size='15')
p1.set_ylabel('Imaginary Axis', size='15')

p1.legend(loc='upper right', prop={'size':10},ncol=1,frameon=True)
p1.set_ylim(0,11)


plt.savefig("vector1.eps",dpi=500)