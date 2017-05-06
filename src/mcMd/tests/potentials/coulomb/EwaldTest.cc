
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/species/Species.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

using namespace McMd;
using namespace Util;

class EwaldTest 
{

public:

   /**
   * Default constructor.
   */
   EwaldTest()
    : nAtom_(0)
   {}

   /**
   * Read parameter file and initialize
   */
   void readParam(std::string paramFileName) 
   {
      // Read parameter file
      std::ifstream in;
      sim.fileMaster().openInputFile(paramFileName, in);
      sim.readParam(in);
      // double diameter;
      // in >> diameter;
      in >> sim.system().boundary();
      in.close();
      std::cout << "Finished reading param file" << std::endl;

      // Allocate and initialize diameters_ array
      int nAtomType = sim.nAtomType();
      diameters_.allocate(nAtomType);
      for (int i = 0; i < nAtomType; ++i) {
         diameters_[i] = 0.2;
      }

      // Allocate and initialize capacities_ array
      int nSpecies = sim.nSpecies();
      capacities_.allocate(nSpecies);
      nAtom_ = 0;
      for (int i = 0; i < nSpecies; ++i) {
         capacities_[i] = sim.species(i).capacity();
         nAtom_ += capacities_[i];
      }

      // Allocate and initialize forces_ array
      forces_.allocate(nAtom_);
      for (int i = 0; i < nAtom_; ++i) {
         forces_[i].zero();
      }

      MdCoulombPotential& coulomb = sim.system().coulombPotential();
      alpha_ = coulomb.get("alpha");
      rSpaceCutoff_ = coulomb.get("rSpaceCutoff");
      kSpaceCutoff_ = coulomb.get("kSpaceCutoff");
   }

   void generateConfig() 
   {
      std::cout << "Begin generating molecules" << std::endl;
      sim.system().generateMolecules(capacities_, diameters_);
      std::cout << "Finished generating molecules" << std::endl;
   };

   void storeForces(DArray<Vector>& forces) 
   {
      UTIL_CHECK(nAtom_ == forces.capacity());
      MdSystem& system = sim.system();
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int iAtom = 0;
      for (int iSpecies = 0; iSpecies < sim.nSpecies(); ++iSpecies) {
         for (system.begin(iSpecies, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               forces[iAtom] = atomIter->force();
               // std::cout << "Atom " << iAtom << "  " << forces[iAtom] << std::endl;
               ++iAtom;
            }
         }
      }
      UTIL_CHECK(iAtom == nAtom_);
   }

   double computeForceError(DArray<Vector> const & refForces) 
   {
      MdSystem& system = sim.system();
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector dF;
      int iAtom = 0;
      double error = 0;
      for (int iSpecies = 0; iSpecies < sim.nSpecies(); ++iSpecies) {
         for (system.begin(iSpecies, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               dF = atomIter->force();
               // std::cout << "atom " << iAtom << " " 
               //           << dF << "  " << refForces[iAtom] << std::endl;
               dF -= refForces[iAtom];
               error += dF.square();
               ++iAtom;
            }
         }
      }
      UTIL_CHECK(iAtom == nAtom_);
      error /= double(nAtom_);
      return sqrt(error);
   }

   void varyAlpha(double alphaMin, double alphaMax, int n)
   {
      MdCoulombPotential& coulomb = sim.system().coulombPotential();
      MdPairPotential& pair = sim.system().pairPotential();

      // Compute well converged system
      coulomb.unsetEnergy();
      pair.unsetEnergy();
      double energyRef = coulomb.energy();
      sim.system().setZeroForces();
      coulomb.addForces();
      pair.addForces();
      storeForces(forces_);

      // Loop over values of alpha
      double dAlpha = (alphaMax - alphaMin)/double(n);
      double alpha, kEnergy, rEnergy, energy, fError;
      for (int i = 0; i <= n; ++i) {
         alpha = alphaMin + dAlpha*i;

         coulomb.set("alpha", alpha);
         coulomb.unsetEnergy();
         pair.unsetEnergy();
         kEnergy = coulomb.kSpaceEnergy();
         rEnergy = coulomb.rSpaceEnergy();
         energy = kEnergy + rEnergy;

         sim.system().setZeroForces();
         coulomb.addForces();
         pair.addForces();
         fError = computeForceError(forces_);

         std::cout << Dbl(coulomb.get("alpha"), 10)
                   << "  " << Dbl(kEnergy, 15)
                   << "  " << Dbl(rEnergy, 15)
                   //<< "  " << Dbl(energy, 20, 13) 
                   << "  " << Dbl(energy - energyRef, 15) 
                   << "  " << Dbl(fError, 15) 
                   << std::endl;
      }

      // Reset to default
      coulomb.set("alpha", alpha_);
   }

   void varyRSpaceCutoff(double rSpaceCutoffMin, double rSpaceCutoffMax, int n)
   {
      MdCoulombPotential& coulomb = sim.system().coulombPotential();
      MdPairPotential& pair = sim.system().pairPotential();

      double dRSpaceCutoff = (rSpaceCutoffMax - rSpaceCutoffMin)/double(n);
      double rSpaceCutoff, kEnergy, rEnergy, energy;
      for (int i = 0; i <= n; ++i) {
         rSpaceCutoff = rSpaceCutoffMin + dRSpaceCutoff*i;
         coulomb.set("rSpaceCutoff", rSpaceCutoff);
         coulomb.unsetEnergy();
         pair.unsetEnergy();
         kEnergy = coulomb.kSpaceEnergy();
         rEnergy = coulomb.rSpaceEnergy();
         energy = kEnergy + rEnergy;
         std::cout << Dbl(coulomb.get("rSpaceCutoff"), 10)
                   << "  " << Dbl(kEnergy, 20)
                   << "  " << Dbl(rEnergy, 20)
                   << "  " << Dbl(energy, 20)
                   << std::endl;
      }

      // Reset to default
      coulomb.set("rSpaceCutoff", rSpaceCutoff_);
   }

   MdSimulation& simulation()
   {  return sim; }

private:

   MdSimulation sim;

   DArray<int> capacities_; 
   DArray<double> diameters_; 
   DArray<Vector> forces_;
   double alpha_;
   double rSpaceCutoff_;
   double kSpaceCutoff_;
   int nAtom_;

};

int main(int argc, char* argv[]) 
{

   EwaldTest test;

   test.readParam("in/param");
   test.generateConfig();
   test.varyAlpha(0.4, 3.0, 26);
   //std::cout << "alpha = " 
   //          << test.simulation().system().coulombPotential().get("alpha")
   //          << std::endl;

}
