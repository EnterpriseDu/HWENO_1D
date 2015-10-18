#from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
import matplotlib.pyplot as plt   
#import matplotlib.animation 
#import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


cut = 0
start = 0
label = 'vortex_GRP9adp'
add = './data/' + label
addx = add + '/xc__'
addy = add + '/yc__'
addr = add + '/rho_'
addi = add + '/rho/'


addlog = './data/' + label + '/log.txt'
file_data = open(addlog)
line = file_data.readline()
line = file_data.readline()
line = file_data.readline()
line = file_data.readline()
line = file_data.readline()
line = file_data.readline()
line = file_data.readline()
file_data.close()

line = line.split(' ')
#N = str.atoi(line[0]) + 1

start = 0
#N = 501
N = 1
for k in range(start, N):
  no = "%04d"%k

  file_name = addx + no + '.txt'
  file_data = open(file_name)
  lines = file_data.readlines()
  file_data.close()
  
  rowx = len(lines)
  datax = []
  for i in range(cut, rowx-cut):
    datax.append([])
    temp = lines[i].split('\t')
    columnx = len(temp)
    for j in range(0, columnx-1):
        datax[i-cut].append(str.atof(temp[j]))
  rowx = len(datax)
  columnx = len(datax[0])

  file_name = addy + no + '.txt'
  file_data = open(file_name)
  lines = file_data.readlines()
  file_data.close()
  rowy = len(lines)
  datay = []
  for i in range(cut, rowy-cut):
    datay.append([])
    temp = lines[i].split('\t')
    columny = len(temp)
    for j in range(0, columny-1):
        datay[i-cut].append(str.atof(temp[j]))

  file_name = addr + no + '.txt'
  file_data = open(file_name)
  lines = file_data.readlines()
  file_data.close()
  row = len(lines)
  datar = []
  for i in range(cut, row-cut):
      datar.append([])
      temp = lines[i].split('\t')#row-1-i
      column = len(temp)
      for j in range(0, column-1):
          datar[i-cut].append(str.atof(temp[j]))
  column = len(datar[0])


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

  file_name = addi + no + '.png'

  cont = plt.contour(datax, datay, datar, 30, cmap = cm.coolwarm)
  fig.colorbar(cont, shrink=0.5, aspect=5)
  plt.savefig(file_name)
  plt.show()
  #del(fig)
  print("===%04d complete==="%k)
