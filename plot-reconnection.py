# Plot reconnection rate

from boutdata import collect
import matplotlib.pyplot as plt
from numpy import pi, argmax

path="merging"
xind=162
zind=320

#path="merging-highres"
#xind=322
#zind=640

#path="merging-lowres"
#xind=82
#zind=160

j = collect("Jpar", path=path, xind=xind, zind=zind)
tau_A = collect("tau_A", path=path)
resistivity = collect("resistivity", path=path)
viscosity = collect("viscosity", path=path)
t_array = collect("t_array", path=path)

mu0 = 4e-7*pi

eta = mu0*resistivity/tau_A # [Ohm m]
print("eta = %e Ohm-m" % eta)
j /= mu0  # Convert to [A/m^2]

rec = -eta * j[:,0,0,0] # -E


# Find maximum reconnection rate
ind = argmax(rec)
print("Maximum reconnection rate: %e V/m at t = %e microseconds" % (rec[ind], t_array[ind]*tau_A*1e6))

plt.plot(t_array, rec*1e-3, label=r"$\eta=%.1e; \mu=%.1e$" % (resistivity, viscosity))
plt.xlabel("Time [t_0 = %.2f microsec]" % (tau_A*1e6))
plt.ylabel(r"Rec. Rate ($\times 10^3$ V/m)")
plt.legend(loc="upper right")

plt.savefig("reconnection-rate.pdf")

plt.show()


# Convert to timescale used in A.Stanier Phys Plasmas 20,122302 (2013)
# Note factor of 2 in normalisation due to B0=0.5 in Stanier, B0=1 here
plt.plot(t_array * tau_A/(0.29e-6), rec*1e-3, label=r"$\eta=%.1e; \mu=%.1e$" % (resistivity*2.0, viscosity*2.0))
plt.xlabel("Time [t_0 = 0.29 microsec]")
plt.ylabel(r"Rec. Rate ($\times 10^3$ V/m)")
plt.legend(loc="upper right")

plt.savefig("reconnection-rate-stanier.pdf")

plt.show()
