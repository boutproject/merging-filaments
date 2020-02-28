
#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

///////////////////////////////////////////////////////
// Add a function coil_greens(Rc, Zc, R, Z) to the
// input expression parser. This calculates the magnetic flux
// due to a single-turn toroidal coil.
 
// Note: C++17 includes elliptic integrals in the standard.
// For now use Cephes library, used by e.g. Scipy
double ellpk(double x); // Complete elliptic integral of the first kind
double ellpe(double x); // Complete elliptic integral of the second kind

/// Calculate poloidal flux at (R,Z) due to a unit current
/// at (Rc,Zc) using Greens function
/// This code was adapted from FreeGS
double coil_greens(double Rc, double Zc, double R, double Z) {
  // Calculate k^2
  double k2 = 4.*R * Rc / ( SQ(R + Rc) + SQ(Z - Zc) );

  // Clip to between 0 and 1 to avoid nans e.g. when coil is on grid point
  if (k2 < 1e-10)
    k2 = 1e-10;
  if (k2 > 1.0 - 1e-10)
    k2 = 1.0 - 1e-10;
  double k = sqrt(k2);

  // Note definition of ellpk, ellpe is K(k^2), E(k^2)
  return 2e-7 * sqrt(R*Rc) * ( (2. - k2)*ellpk(k2) - 2.*ellpe(k2) ) / k;
}

class CoilGenerator : public FieldGenerator {
public:
  CoilGenerator() = default;
  
  CoilGenerator(BoutReal Rc, BoutReal Zc, FieldGeneratorPtr R, FieldGeneratorPtr Z)
    : Rc(Rc), Zc(Zc), Rgen(R), Zgen(Z) {}
  
  BoutReal generate(BoutReal x, BoutReal y, BoutReal z, BoutReal t) override {
    // Calculate R,Z location of this (x,y,z) point
    BoutReal R = Rgen->generate(x,y,z,t);
    BoutReal Z = Zgen->generate(x,y,z,t);
    
    return coil_greens(Rc, Zc, R, Z);
  }
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() != 4) {
      throw BoutException("coil_greens expects 4 arguments (Rc, Zc, R, Z)");
    }
    auto argsit = args.begin();
    // Coil positions are constants
    BoutReal Rc_new = (*argsit++)->generate(0,0,0,0);
    BoutReal Zc_new = (*argsit++)->generate(0,0,0,0);
    // Evaluation location can be a function of x,y,z,t
    FieldGeneratorPtr Rgen_new = *argsit++;
    FieldGeneratorPtr Zgen_new = *argsit;
    return std::make_shared<CoilGenerator>(Rc_new, Zc_new, Rgen_new, Zgen_new);
  }
private:
  BoutReal Rc, Zc;  // Coil location, fixed
  FieldGeneratorPtr Rgen, Zgen; // Location, function of x,y,z,t
};

///////////////////////////////////////////////////////

class MergingFlux : public PhysicsModel {
protected:
  int init(bool restarting) {

    // Add a function which can be used in input expressions
    // This calculates the Greens function for a coil
    FieldFactory::get()->addGenerator("coil_greens", std::make_shared<CoilGenerator>());
    
    // Get the magnetic field
    Coordinates *coord = mesh->getCoordinates();
    B0 = coord->Bxy;
    R = sqrt(coord->g_22);
    
    // Read options
    Options& opt = Options::root()["model"];
    resistivity = opt["resistivity"].withDefault(0.0);
    viscosity = opt["viscosity"].withDefault(0.0);
    density = opt["density"].withDefault(5e18);
    Te = opt["Te"].doc("Reference temperature [eV]").withDefault(10.);
    AA = opt["AA"].doc("Atomic mass number").withDefault(2.0);

    Bv = opt["Bv"].doc("External vertical field [T]").withDefault(0.0);
    Psiext = opt["Psiext"].doc("External magnetic flux [Wb]").withDefault(Field3D(0.0));

    SAVE_ONCE(Psiext);
    
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
  int rhs(BoutReal UNUSED(t)) {
    
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
    psi = Apar * R + 0.5*Bv*SQ(R) + Psiext;
    
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
  Field3D Psiext; // External magnetic flux
};

BOUTMAIN(MergingFlux);
