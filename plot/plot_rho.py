#from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
import matplotlib.pyplot as plt
from matplotlib import cm


cut = 0
fname = '0000'
label = 'test'
add = './data/' + label
#addx = add + '/xc__'
#addy = add + '/yc__'
addr = add + '/rho_'
#addx = addx + fname
#addy = addy + fname
addr = addr + fname
#addx = addx + '.txt'
#addy = addy + '.txt'
addr = addr + '.txt'







file = open(addr)
data_str = file.readline()
file.close()
datar = []
temp = data_str.split('\t')#row-1-i
column = len(temp)
for j in range(0, column-1):
    datar.append(str.atof(temp[j]))
column = len(datar)

xl = 0.0
xr = 100.0

h_x = (xr - xl)/column;
x1 = xl + 0.5 * h_x;
xm = xr - 0.5 * h_x;
x = (np.linspace(x1,xm,column))


"""
"""
fig = plt.figure()

plt.plot(x, datar, '*')

plt.show()
#plt.savefig('./data_out/test.png') 
