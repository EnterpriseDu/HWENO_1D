#from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
import matplotlib.pyplot as plt   
#import matplotlib.animation 
#import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



fname = '0000'
label = 'vortex_GRP8adp'
add = './data/' + label
addx = add + '/x___'
addy = add + '/y___'
addx = addx + fname
addy = addy + fname
addx = addx + '.txt'
addy = addy + '.txt'


myfilter = 1
file = open(addx)
linex = file.readlines()
file.close()
rowx = len(linex)
datax = []
for i in range(0, rowx):
  if(i%myfilter == 0):
    datax.append([])
    temp = linex[i].split('\t')
    columnx = len(temp)
    for j in range(0, columnx-1):
      if(j%myfilter == 0):
        datax[i/myfilter].append(str.atof(temp[j]))
rowx = len(datax)
columnx = len(datax[0])
print("row=%d, column=%d"%(rowx, columnx))


file = open(addy)
liney = file.readlines()
file.close()
rowy = len(liney)
datay = []
for i in range(0, rowy):
  if(i%myfilter == 0):
    datay.append([])
    temp = liney[i].split('\t')
    columny = len(temp)
    for j in range(0, columny-1):
      if(j%myfilter == 0):
        datay[i/myfilter].append(str.atof(temp[j]))

print("row=%d, column=%d"%(rowx, columnx))

datax2 = []
datay2 = []
for i in range(0, columnx):
  datax2.append([])
  datay2.append([])

for i in range(0, columnx):
  for j in range(0, rowx):
    datax2[i].append(datax[j][i])
    datay2[i].append(datay[j][i])


fig1 = plt.figure(figsize=(8,8))

xmajorLocator   = MultipleLocator(20.0) 
xmajorFormatter = FormatStrFormatter('%d')
#xminorLocator   = MultipleLocator(5)
ymajorLocator   = MultipleLocator(20.0)
ymajorFormatter = FormatStrFormatter('%d')
#yminorLocator   = MultipleLocator(0.1)
ax = plt.subplot(111)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_major_formatter(ymajorFormatter)
ax.xaxis.grid(True, which='major')  
ax.yaxis.grid(True, which='major')

for i in range(0, rowx):
      plt.plot(datax[i], datay[i], color='black')
for i in range(0, columnx):
      plt.plot(datax2[i], datay2[i], color='black')

plt.show()
