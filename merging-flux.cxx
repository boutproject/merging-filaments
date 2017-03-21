
#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <invert_laplace.hxx>

class MergingFlux : public PhysicsModel {
protected:
  int init(bool restarting) {
    
    // Get the magnetic field
    Coordinates *coord = mesh->coordinates();
    B0 = coord->Bxy;
    R = sqrt(coord->g_22);
    
    // Read options
    Options *opt = Options::getRoot()->getSection("model");
    OPTION(opt, resistivity, 0.0);
    OPTION(opt, viscosity, 0.0);
    OPTION(opt, density, 5e18);
    OPTION(opt, Te, 10.); // Reference temperature [eV]
    OPTION(opt, AA, 2.0); // Atomic mass number

    OPTION(opt, Bv, 0.0); // External vertical field [T]
    
    BoutReal mass_density = density * AA*SI::Mp; // kg/m^3

    // Calculate normalisations
    tau_A = sqrt(SI::mu0 * mass_density); // timescale [s]
    
    SAVE_ONCE2(tau_A, density); // Save to output
    SAVE_ONCE2(resistivity, viscosity); 

    output.write("\tDensity = %e m^-3, tau_A = %e seconds\n", density, tau_A);

    // Collisional calculation
    
    BoutReal Coulomb = 6.6 - 0.5*log(density/1e20) + 1.5*log(Te); // Coulomb logarithm
    output.write("\tCoulomb logarithm = %e\n", Coulomb);
    
    // Electron collision time [seconds]
    BoutReal tau_e = 1. / (2.91e-6 * (density / 1e6) * Coulomb * pow(Te, -3./2));
    // ion-ion collision time [seconds]
    BoutReal tau_i = sqrt(AA) / (4.80e-8 * (density / 1e6) * Coulomb * pow(Te, -3./2));
    
    output.write("\tCollision times [sec]: tau_e = %e, tau_i = %e\n", tau_e, tau_i);
    
    // Parallel conductivity
    BoutReal sigma_par = 1.96*density*SQ(SI::qe)*tau_e/SI::Me;
    
    output.write("\tBraginskii resistivity: %e [Ohm m]\n", 1./sigma_par);
    output.write("\tNormalised Braginskii: %e\n", tau_A/(SI::mu0*sigma_par));
    output.write("\tUsing resistivity: %e\n", resistivity);
    
    
    // Perpendicular viscosity
    BoutReal Bmag = max(B0, true); // max over domain
    BoutReal Omega_i = SI::qe*Bmag/(AA*SI::Mp);
    BoutReal eta_perp = 0.3 * density*SI::qe*Te/ ( SQ(Omega_i)*tau_i );
    
    // Perpendicular gyro-viscosity
    BoutReal eta_gyro = 0.5*density*SI::qe*Te/Omega_i;
    
    output.write("\tViscosity: %e [Pa s], Gyro-viscosity: %e [Pa s]\n", eta_perp, eta_gyro);
    output.write("\tKinematic viscosity: %e [m^2/s], Gyro-viscosity: %e [m^2/s]\n", eta_perp/mass_density, eta_gyro/mass_density);
    output.write("\tNormalised kin. viscosity: %e, gyro: %e\n", tau_A*eta_perp/mass_density, tau_A*eta_gyro/mass_density);
    output.write("\tUsing viscosity: %e\n", viscosity);

    // Solve for electromagnetic potential and vorticity
    SOLVE_FOR2(Apar, omega);
    
    // Read parallel current density
    mesh->get(Jpar, "Jpar");

    // Create Laplacian inversion objects for potentials
    phiSolver = Laplacian::create(Options::getRoot()->getSection("phisolver"));
    psiSolver = Laplacian::create(Options::getRoot()->getSection("psisolver"));

    if(!restarting) {
      // Invert Jpar to get vector potential
      
      Apar = -psiSolver->solve(Jpar);
    }
    phi = 0.0; // Initial value
    
    // Additional outputs
    SAVE_REPEAT3(phi, Jpar, psi);
    
    return 0;
  }
  int rhs(BoutReal t) {
    
    mesh->communicate(Apar, omega);
    
    // Apply Dirichlet boundary conditions in z
    for(int i=0;i<mesh->LocalNx;i++) {
      for(int j=0;j<mesh->LocalNy;j++) {
        Apar(i,j,0) = -Apar(i,j,1);
        Apar(i,j,mesh->LocalNz-1) = -Apar(i,j,mesh->LocalNz-2);
        
        omega(i,j,0) = -omega(i,j,1);
        omega(i,j,mesh->LocalNz-1) = -omega(i,j,mesh->LocalNz-2);
      }
    }

    // Get J from Psi
    Jpar = -Delp2(Apar);

    // Get phi from vorticity
    phi = phiSolver->solve(omega*SQ(B0));
    mesh->communicate(Jpar, phi);

    Jpar.applyBoundary("neumann");
    for(int i=0;i<mesh->LocalNx;i++) {
      for(int j=0;j<mesh->LocalNy;j++) {
        Jpar(i,j,0) = Jpar(i,j,1);
        Jpar(i,j,mesh->LocalNz-1) = Jpar(i,j,mesh->LocalNz-2);
      }
    }
    
    // Poloidal flux, including external field
    psi = Apar * R + 0.5*Bv*SQ(R);
    
    // Vorticity
    ddt(omega) = 
      - bracket(psi, Jpar, BRACKET_ARAKAWA)/R // b dot Grad(Jpar)
      - bracket(phi, omega, BRACKET_ARAKAWA)  // ExB advection
      + viscosity*Delp2(omega)  // Viscosity
      ;
    
    // Vector potential
    ddt(Apar) = 
      bracket(psi, phi, BRACKET_ARAKAWA)/R // b dot Grad(phi)
      - resistivity * Jpar // Resistivity
      ;
    
    return 0;
  }

private:
  // Evolving variables
  Field3D Apar, omega;  // Electromagnetic potential, vorticity
  
  Field3D Jpar;   // Parallel current density
  Field3D phi;    // Electrostatic potential
  Field3D psi;    // Poloidal flux

  Field2D B0;  // Magnetic field [T]

  Laplacian *phiSolver; // Solver for potential phi from vorticity
  Laplacian *psiSolver; // Solver for psi from current Jpar

  BoutReal resistivity;
  BoutReal viscosity;
  
  BoutReal tau_A; // Normalisation timescale [s]
  BoutReal density; // Normalisation density [m^-3]
  BoutReal Te; // Temperature for collisions calculation [eV]
  BoutReal AA; // Atomic mass number

  Field2D R; // Major radius
  BoutReal Bv; // Vertical magnetic field
};

BOUTMAIN(MergingFlux);
