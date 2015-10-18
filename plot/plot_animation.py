from mpl_toolkits.axes_grid.axislines import SubplotZero
import string as str
import numpy as np   
import matplotlib.pyplot as plt   
import matplotlib.animation 
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



xl = 0.0
xr = 100.0
yb = 0.0
yt = 100.0
sol = './data_out/rho_200.txt'
log = './data_out/log_200.txt'



file = open(log)
lines = file.readlines()
str_N = lines[9][8:len(lines[9])-1]
N = str.atoi(str_N)
file.close()
lines = []
print("%d" %N)

file = open(sol)
lines = file.readlines()
file.close()
row = len(lines)
idx = 0;
n = row/(N+1)
print("%d" %row)
print("%d" %n)



data = []
for k in range(0, N+1):
  data.append([])
  for i in range(0, n):
    data[k].append([])
    temp = lines[k*n+i].split('\t')
    column = len(temp)
    for j in range(0, column-1):
      data[k][i].append(str.atof(temp[j]))

column = len(data[0][0])



#for k in range(0, N):
#  for i in range(0, n):
#    for j in range(0, column):
#      print(data[k][i][j]),
#      print(' '),
#    print('\n'),
#  print('\n'),


h_x = (xr - xl)/column;
h_y = (yt - yb) / n;
x1 = xl + 0.5 * h_x;
xm = xr - 0.5 * h_x;
y1 = yb + 0.5 * h_y;
yn = yt - 0.5 * h_y;
x = (np.linspace(x1,xm,column))
y = (np.linspace(y1,yn,n))
X,Y = np.meshgrid(x,y)


"""
#=====================ANIMATION 3D CONTOUR/SURFACE=========================
fig = plt.figure()
ax = axes3d.Axes3D(fig)


#wframe = ax.plot_wireframe(X, Y, data[0], cmap=cm.coolwarm)
wframe = ax.contour(X, Y, data[i], extend3d = True, cmap=cm.coolwarm)
#wframe = ax.scatter(X, Y, data[0], c='b', marker='o')
#wframe = ax.plot_surface(X, Y, data[0], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_zlim(-1.2,1.2)
fig.colorbar(wframe, shrink=0.5, aspect=5)



def update(i, ax, fig):
    ax.cla()
    #wframe = ax.plot_wireframe(X, Y, data[i], cmap=cm.coolwarm)
    wframe = ax.contour(X, Y, data[i], extend3d = True, cmap=cm.coolwarm)
    #wframe = ax.scatter(X, Y, data[i], c='b', marker='o')
    #wframe = ax.plot_surface(X, Y, data[i], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(-1.2,1.2)
    return wframe,

ani = matplotlib.animation.FuncAnimation(fig, update, frames=N, fargs=(ax, fig), interval=100)
"""



#========================ANIMATION 2D CONTOUR==============================
fig = plt.figure()
"""
xmajorLocator   = MultipleLocator(0.5) 
xmajorFormatter = FormatStrFormatter('%1.1f')
#xminorLocator   = MultipleLocator(5)
ymajorLocator   = MultipleLocator(0.5)
ymajorFormatter = FormatStrFormatter('%1.1f')
#yminorLocator   = MultipleLocator(0.1)
"""
ax = plt.subplot(111)
"""
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_major_formatter(ymajorFormatter)
ax.xaxis.grid(True, which='major')  
ax.yaxis.grid(True, which='major')
"""
cont = plt.contour(X, Y, data[0], 30, cmap = cm.coolwarm)
fig.colorbar(cont, shrink=0.5, aspect=5)


# animation function
def animate(i): 
    ax.cla()
    #ax.xaxis.set_major_locator(xmajorLocator)
    #ax.xaxis.set_major_formatter(xmajorFormatter)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_major_formatter(ymajorFormatter)
    #ax.xaxis.grid(True, which='major')  
    #ax.yaxis.grid(True, which='major')
    cont = plt.contour(X, Y, data[i], 30, cmap = cm.coolwarm)
    return cont  

anim = matplotlib.animation.FuncAnimation(fig, animate, frames=N, interval = 10)



#plt.contour(X, Y, data[15], 9, cmap = cm.coolwarm)

plt.show()

"""
print(N)
print(column)
print(h_x)
print(h_y)
print(yb)
print(yt)
print(xl)
print(xr)
"""
