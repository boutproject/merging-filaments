
from boutdata import collect
from numpy import argmin, abs, pi, sqrt, arange, interp, transpose, linspace, concatenate, amin, amax
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker

xlim = None #[0.3, 1.5] # Radial limits [m]
ylim = None #[-0.6, +0.6] # Vertical limits [m]

path="toroidal"

t_array = collect("t_array", path=path)
tau_A = collect("tau_A", path=path)

dx = collect("dx", path=path)[0,0]
g_11 = collect("g_11", path=path)[0,0]
dz = collect("dz", path=path)
g_33 = collect("g_33", path=path)[0,0]

dx = dx * sqrt(g_11) # Length in x
dz = dz * sqrt(g_33) # Length in z

time = t_array*tau_A*1e6 # in microseconds
nt = len(time) # Number of time points

mu0 = 4e-7*pi

jmin = None
jmax = None

psimin = None
psimax = None

tinds = range(nt)

for tind in tinds:
    # Read current density
    j = collect("Jpar", tind=tind, path=path)[0,:,0,:] / mu0
    
    # Read psi
    psi = collect("psi", tind=tind, path=path)[0,:,0,:]

    if jmin is None:
        jmin = amin(j)
        jmax = amax(j)
        psimin = amin(psi)
        psimax = amax(psi)
    else:
        jmin = min(jmin, amin(j))
        jmax = max(jmax, amax(j))
        psimin = min(psimin, amin(psi))
        psimax = max(psimax, amax(psi))


nx,nz = j.shape

x = (arange(nx)-1.5)*dx + 0.2
z = (arange(nz)+0.5)*dz - 2.2

# Now have ranges. Define levels for plotting

alevels = linspace(psimin, psimax, 20)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
cmap = plt.cm.get_cmap("RdBu_r")

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

def fmt(x, pos):
    """
    Format number as LaTeX in scientific notation
    Taken from: http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
    """
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

for tind in tinds:
    # Read current density
    j = collect("Jpar", tind=tind, path=path)[0,:,0,:] / mu0
    
    # Read psi
    psi = collect("psi", tind=tind, path=path)[0,:,0,:]
    
    ax.clear()
    c = ax.contourf(x,z,transpose(j), 401,
                    norm=MidpointNormalize(vmin=jmin,
                                           vmax=jmax,
                                           midpoint=0.),
                    cmap=cmap)

    if tind == tinds[0]:
        # Only add a color bar once
        fig.colorbar(c, ax=ax, extend='both', format=ticker.FuncFormatter(fmt))

    ax.contour(x, z, transpose(psi), alevels, colors='k')

    ax.set_aspect('equal')

    plt.xlabel("R [m]")
    plt.ylabel("Z [m]")
    plt.title(r"$J_{||}$ [A/m$^2$] at Time: $%.2f\mu$s" % time[tind])
    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.savefig("toroidal_jpar_t%03d.jpg" % (tind), bbox_inches='tight', pad_inches = 0)

plt.show()
