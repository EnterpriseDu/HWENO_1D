import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation

# This example uses subclassing, but there is no reason that the proper function
# couldn't be set up and then use FuncAnimation. The code is long, but not
# really complex. The length is due solely to the fact that there are a total
# of 9 lines that need to be changed for the animation as well as 3 subplots
# that need initial set up.
class SubplotAnimation(animation.TimedAnimation):
    def __init__(self, dataRHO, dataU, dataP, l, r, column, row):
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 4)

        #self.t = np.linspace(0, 80, 400)
        #self.x = np.cos(2 * np.pi * self.t / 10.)
        #self.y = np.sin(2 * np.pi * self.t / 10.)
        #self.z = 10 * self.t

        self.x = np.linspace(l,r,column)
        self.t = np.linspace(0,100,row)
        self.RHO = dataRHO
        self.U = dataU
        self.P = dataP

        ax1.set_xlabel('x')
        ax1.set_ylabel('RHO')
        self.line1 = Line2D([], [])
        #self.line1a = Line2D([], [], color='red', linewidth=2)
        #self.line1e = Line2D([], [], color='red', marker='o', markeredgecolor='r')
        ax1.add_line(self.line1)
        #ax1.add_line(self.line1a)
        #ax1.add_line(self.line1e)
        #ax1.set_xlim(-1, 1)
        #ax1.set_ylim(-2, 2)
        #ax1.set_aspect('equal', 'datalim')

        ax2.set_xlabel('x')
        ax2.set_ylabel('U')
        self.line2 = Line2D([], [])
        #self.line2a = Line2D([], [], color='red', linewidth=2)
        #self.line2e = Line2D([], [], color='red', marker='o', markeredgecolor='r')
        ax2.add_line(self.line2)
        #ax2.add_line(self.line2a)
        #ax2.add_line(self.line2e)
        #ax2.set_xlim(-1, 1)
        #ax2.set_ylim(0, 800)

        ax3.set_xlabel('x')
        ax3.set_ylabel('P')
        self.line3 = Line2D([], [])
        #self.line3a = Line2D([], [], color='red', linewidth=2)
        #self.line3e = Line2D([], [], color='red', marker='o', markeredgecolor='r')
        ax3.add_line(self.line3)
        #ax3.add_line(self.line3a)
        #ax3.add_line(self.line3e)
        #ax3.set_xlim(-1, 1)
        #ax3.set_ylim(0, 800)

        animation.TimedAnimation.__init__(self, fig, interval=50, blit=True)

    def _draw_frame(self, framedata):
        i = framedata
        head = i - 1
        head_len = 10
        head_slice = (self.t > self.t[i] - 1.0) & (self.t < self.t[i])

        self.line1.set_data(self.x, self.RHO[i])
        #self.line1a.set_data(self.x[head_slice], self.y[head_slice])
        #self.line1e.set_data(self.x[head], self.y[head])

        self.line2.set_data(self.x, self.U[i])
        #self.line2a.set_data(self.y[head_slice], self.z[head_slice])
        #self.line2e.set_data(self.y[head], self.z[head])

        self.line3.set_data(self.x, self.P[i])
        #self.line3a.set_data(self.x[head_slice], self.z[head_slice])
        #self.line3e.set_data(self.x[head], self.z[head])

        self._drawn_artists = [self.line1, self.line2, self.line3]

    def new_frame_seq(self):
        return iter(range(self.t.size))

    def _init_draw(self):
        lines =  [self.line1, self.line2, self.line3]
        for l in lines:
            l.set_data([], [])

ani = SubplotAnimation()
#ani.save('test_sub.mp4')
plt.show()
