
from boutdata import collect
from numpy import argmin, abs, pi, sqrt, arange, interp, transpose
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Time in microseconds where a plot is desired
times = [7.54, 10.54, 18.04]

xlim = [0.3, 1.5] # Radial limits [m]
ylim = [-0.6, +0.6] # Vertical limits [m]

path="merging"

t_array = collect("t_array", path=path)
tau_A = collect("tau_A", path=path)

dx = collect("dx", path=path)[0,0]
g_11 = collect("g_11", path=path)[0,0]
dz = collect("dz", path=path)
g_33 = collect("g_33", path=path)[0,0]

dx = dx * sqrt(g_11) # Length in x
dz = dz * sqrt(g_33) # Length in z

time = t_array*tau_A*1e6 # in microseconds

mu0 = 4e-7*pi

class MidpointNormalize(colors.Normalize):
    """
    Taken from http://matplotlib.org/users/colormapnorms.html
    From Joe Kington: This one gives two different linear ramps:
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

for t in times:
    # Find nearest time index
    tind = argmin(abs(time - t))
    
    # Read current density
    j = collect("Jpar", tind=tind, path=path)[0,:,0,:] / mu0
    
    # Read Apar
    Apar = collect("Apar", tind=tind, path=path)[0,:,0,:]

    nx,nz = j.shape
    
    x = (arange(nx)+0.5)*dx
    z = (arange(nz)+0.5)*dz - 2.2
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    pcm = ax.pcolormesh(x, z, transpose(j),
                           norm=MidpointNormalize(midpoint=0.),
                           cmap='RdBu_r')
    fig.colorbar(pcm, ax=ax, extend='both')

    ax.contour(x, z, transpose(Apar), 20, colors='k')

    ax.set_aspect('equal')

    plt.xlabel("R [m]")
    plt.ylabel("Z [m]")
    plt.title(r"$J_{||}$ [A/m$^2$] at Time: $%.2f\mu$s" % time[tind])
    plt.xlim(xlim)
    plt.ylim(ylim)

    #plt.savefig("jpar_t%02d.pdf" % (int(time[tind])))
    plt.savefig("jpar_t%02d.png" % (int(time[tind])))
    plt.show()
