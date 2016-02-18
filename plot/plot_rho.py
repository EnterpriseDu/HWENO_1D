#from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
import matplotlib.pyplot as plt
from matplotlib import cm



fname = '0000'
label = 'wood'
scheme = 'G4H5'
version = 'dev-1.1'
m = '800'
switches = 'FWDY'
label = label + '_' + scheme + '_' + version + '_' + m + '_' + switches
#label = 'DM_GRP-H-fixed_limiter-test_240x60_002010'
add = '../SOLUTION/' + label
addr = add + '/rho_'
addr = addr + fname
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
