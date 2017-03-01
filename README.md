Merging current filaments
=========================

2D simulation in X-Z of merging current filaments
in a rectangular domain with ideal walls. 

Evolves flute-reduced, zero-beta model for vorticity U
and magnetic flux Psi, with auxilliary variables Jpar
(parallel current density) and phi (electrostatic potential).

    dU/dt = [Psi, Jpar] - [phi, U] + nu*Delp2(U)
    dPsi/dt = [Psi, phi] + eta * Jpar

Test cases are

*single*  A single current filament at the centre of the domain
          with eta = 1e-5

*single-no-resist*  A single filament with eta=0

*merging*  Two current filaments with eta=1e-5

*merging-no-resist*  Two current filaments with eta = 0

Run all cases, then the python script "makeplots.py"

    $ make
    $ mpirun -np 8 ./merging-flux -d single
    $ mpirun -np 8 ./merging-flux -d single-no-resist
    $ mpirun -np 8 ./merging-flux -d merging
    $ mpirun -np 8 ./merging-flux -d merging-no-resist
    $ python makeplots.py

This should produce a number of PDF figures

MAST simulations
----------------

Typical parameters prior to merging:

density n_0 = 5e18 m^-3
Magnetic field B_0 = 0.5 T
Te = Ti = 10eV



