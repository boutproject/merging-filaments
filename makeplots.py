
import matplotlib.pyplot as plt

from boutdata import collect

from numpy import amax, abs, linspace, transpose

# Single filament results

psi1 = collect("psi", path="single", xind=40, zind=80)[:,0,0,0]
psi2 = collect("psi", path="single-no-resist", xind=40, zind=80)[:,0,0,0]
t = collect("t_array", path="single")

plt.plot(t, psi1, label=r"$\eta = 10^{-5}$")
plt.plot(t, psi2, label=r"$\eta = 0.0$")
plt.xlabel("Time [normalised]")
plt.ylabel(r"Flux $\psi$ at filament centre")
plt.legend()
plt.savefig("single_flux.pdf")
plt.show()

# Shape change
j1 = collect("Jpar", path="single", zind=80,tind=0)[0,:,0,0]
j2 = collect("Jpar", path="single", zind=80,tind=1)[0,:,0,0]
j3 = collect("Jpar", path="single", zind=80,tind=100)[0,:,0,0]

plt.plot(j1, label=r"$t=0$")
plt.plot(j2, label=r"$t=5$")
plt.plot(j3, label=r"$t=500$")

plt.xlabel("X grid index")
plt.ylabel(r"$J_{||}$")
plt.legend()

plt.savefig("single_jpar_x.pdf")
plt.show()

j1 = collect("Jpar", path="single", xind=40,tind=0)[0,0,0,:]
j2 = collect("Jpar", path="single", xind=40,tind=1)[0,0,0,:]
j3 = collect("Jpar", path="single", xind=40,tind=100)[0,0,0,:]

plt.plot(j1, label=r"$t=0$")
plt.plot(j2, label=r"$t=5$")
plt.plot(j3, label=r"$t=500$")

plt.xlabel("Z grid index")
plt.ylabel(r"$J_{||}$")
plt.legend()

plt.savefig("single_jpar_z.pdf")
plt.show()

### Merging

psi1 = collect("psi", path="merging", xind=40, zind=80)[:,0,0,0]
psi2 = collect("psi", path="merging-no-resist", xind=40, zind=80)[:,0,0,0]
t = collect("t_array", path="single")

plt.plot(t, psi1, label=r"$\eta = 10^{-5}$")
plt.plot(t, psi2, label=r"$\eta = 0.0$")
plt.xlabel("Time [normalised]")
plt.ylabel(r"Flux $\psi$ at domain centre")
plt.legend()
plt.savefig("merge_flux.pdf")
plt.show()

def contourJ(data, nlevels=51):
    maxabs = amax(abs(data))
    levels = linspace(-maxabs, maxabs, num=nlevels)
    
    cmap = plt.cm.get_cmap("bwr")
    plt.contourf(transpose(data), levels, cmap=cmap)
    plt.colorbar()
    plt.xlabel("X grid index")
    plt.ylabel("Z grid index")

for tind in [0, 100]:
    j = collect("jpar", path="merging-no-resist", tind=tind)

    contourJ(j[0,:,0,:])
    
    plt.title(r"$J_{||}$ at t = %03d. $\eta = 0$" % (t[tind],))

    plt.savefig("merge-nores-j-t%03d.pdf" % (tind,))
    plt.show()

for tind in [0, 20, 40, 100]:
    j = collect("jpar", path="merging", tind=tind)

    contourJ(j[0,:,0,:])
    
    plt.title(r"$J_{||}$ at t = %03d. $\eta = 10^{-5}$" % (t[tind],))

    plt.savefig("merge-j-t%03d.pdf" % (tind,))
    plt.show()
