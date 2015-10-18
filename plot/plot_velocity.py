#from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
#import matplotlib.pyplot as plt   
#import matplotlib.animation 
#import mpl_toolkits.mplot3d.axes3d as p3
#from matplotlib import cm
#from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from pylab import *


fname = '0606'
label = 'inter_shock'
add = './data_out/' + label
addx = add + '/xc__'
addy = add + '/yc__'
addu = add + '/u___'
addv = add + '/v___'
addx = addx + fname
addy = addy + fname
addu = addu + fname
addv = addv + fname
addx = addx + '.txt'
addy = addy + '.txt'
addu = addu + '.txt'
addv = addv + '.txt'

step = 4

file = open(addx)
lines = file.readlines()
file.close()
rowx = len(lines)
datax = []
for i in range(0, rowx):
  if(not(i%step)):
    datax.append([])
    temp = lines[i].split('\t')
    columnx = len(temp)
    for j in range(0, columnx-1):
      if(not(j%step)):
        datax[i/step].append(str.atof(temp[j]))
rowx = len(datax)
columnx = len(datax[0])
print("row=%d, column=%d"%(rowx, columnx))


file = open(addy)
lines = file.readlines()
file.close()
rowy = len(lines)
datay = []
for i in range(0, rowy):
  if(not(i%step)):
    datay.append([])
    temp = lines[i].split('\t')
    columny = len(temp)
    for j in range(0, columny-1):
      if(not(j%step)):
        datay[i/step].append(str.atof(temp[j]))


file = open(addu)
lines = file.readlines()
file.close()
row = len(lines)
datau = []
for i in range(0, row):
  if(not(i%step)):
    datau.append([])
    temp = lines[i].split('\t')#row-1-i
    column = len(temp)
    for j in range(0, column-1):
      if(not(j%step)):
        datau[i/step].append(str.atof(temp[j]))


file = open(addv)
lines = file.readlines()
file.close()
row = len(lines)
datav = []
for i in range(0, row):
  if(not(i%step)):
    datav.append([])
    temp = lines[i].split('\t')#row-1-i
    column = len(temp)
    for j in range(0, column-1):
      if(not(j%step)):
        datav[i/step].append(str.atof(temp[j]))


row = len(datau)
datar = []
for i in range(0, row):
    datar.append([])
    column = len(datau[i])
    for j in range(0, column-1):
        datar[i].append(np.sqrt(datau[i][j]**2+datav[i][j]**2))
column = len(datav[0])


fig = plt.figure()

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

quiv = quiver(datax,datay,datau,datav, datar, cmap = cm.coolwarm)
fig.colorbar(quiv, shrink=0.5, aspect=5)

plt.show()

