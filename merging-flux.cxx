
#include <bout/physicsmodel.hxx>
#include <invert_laplace.hxx>

class MergingFlux : public PhysicsModel {
protected:
  int init(bool restarting) {
    // Read options
    Options *opt = Options::getRoot()->getSection("model");
    OPTION(opt, resistivity, 0.0);
    OPTION(opt, viscosity, 0.0);
    OPTION(opt, density, 5e18);
    
    // Calculate normalisations
    tau_alfven = sqrt(4e-7*PI * density);
    
    SOLVE_FOR2(Psi, U);
    
    // Read parallel current density
    mesh->get(Jpar, "Jpar");

    // Create Laplacian inversion objects for potentials
    phiSolver = Laplacian::create(Options::getRoot()->getSection("phisolver"));
    psiSolver = Laplacian::create(Options::getRoot()->getSection("psisolver"));

    if(!restarting) {
      // Invert Jpar to get vector potential
      
      Psi = psiSolver->solve(Jpar);
    }
    phi = 0.0; // Initial value
    
    // Additional outputs
    SAVE_REPEAT2(phi, Jpar);

    return 0;
  }
  int rhs(BoutReal t) {
    
    Coordinates *coord = mesh->coordinates();
    mesh->communicate(Psi, U);

    // Get J from Psi
    Jpar = Delp2(Psi);

    // Get phi from vorticity
    phi = phiSolver->solve(U);
    mesh->communicate(Jpar, phi);

    // Vorticity
    ddt(U) = 
      SQ(coord->Bxy) * bracket(Psi, Jpar, BRACKET_ARAKAWA) // b dot Grad(Jpar)
      - bracket(phi, U, BRACKET_ARAKAWA)  // ExB advection
      + viscosity*Delp2(U)  // Viscosity
      ;
    
    // Vector potential
    ddt(Psi) = 
      bracket(Psi, phi, BRACKET_ARAKAWA) // b dot Grad(phi)
      + resistivity * Jpar // Resistivity
      ;
    
    return 0;
  }

private:
  // Evolving variables
  Field3D Psi, U;  // Electromagnetic potential, vorticity
  
  Field3D Jpar;   // Parallel current density
  Field3D phi;    // Electrostatic potential

  Laplacian *phiSolver; // Solver for potential phi from vorticity
  Laplacian *psiSolver; // Solver for psi from current Jpar

  BoutReal resistivity;
  BoutReal viscosity;
  BoutReal density;
};

BOUTMAIN(MergingFlux);
