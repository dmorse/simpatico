#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
#include <mcMd/potentials/coulomb/MdSpmePotential.h>
#include <mcMd/potentials/coulomb/MdEwaldPotential.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <simp/species/Species.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

using namespace McMd;
using namespace Util;

class SpmeTest 
{

public:

   /**
   * Default constructor.
   */
   SpmeTest()
    : ewald(sim.system()),
      nAtom_(0)
   {
   }

   /**
   * Read parameter file and initialize
   */
   void readParam(std::string paramFileName) 
   {
      // Read parameter file
      std::ifstream in;
      sim.fileMaster().openInputFile(paramFileName, in);
      sim.readParam(in);
      ewald.readParam(in);
      ewald.rSpaceAccumulator().setPairPotential(sim.system().pairPotential());
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

      alpha_ = sim.system().coulombPotential().get("alpha");
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

   void compareKSpace()
   {
      MdCoulombPotential& spme = sim.system().coulombPotential();

      // Compute well converged system
      ewald.unsetEnergy();
      ewald.computeEnergy();
      double kEnergyRef = ewald.kSpaceEnergy();
      kEnergyRef /= double(nAtom_);

      sim.system().setZeroForces();
      ewald.addForces();
      storeForces(forces_);

      spme.unsetEnergy();
      spme.makeWaves();
      spme.computeEnergy();
      double kEnergy = spme.kSpaceEnergy();
      kEnergy /= double(nAtom_);

      std::cout << "kEnergy= " << Dbl(kEnergy, 15) 
                << std::endl;
      std::cout << "eError = "<< Dbl(kEnergy - kEnergyRef, 15) 
                << std::endl;

      sim.system().setZeroForces();
      spme.addForces();
      double fError = computeForceError(forces_);
      std::cout << "fError = " << Dbl(fError, 15) << std::endl;
   }

   void varyAlpha(double alphaMin, double alphaMax, int n)
   {
      MdCoulombPotential& spme = sim.system().coulombPotential();
      MdPairPotential& pair = sim.system().pairPotential();
      double volume = sim.system().boundary().volume();
      double nAtom = double(nAtom_);

      // Set Ewald and SPME to same parameters
      ewald.set("alpha", 1.4);
      spme.set("alpha", 1.4);

      // Compute reference energy, well converged system
      ewald.unsetEnergy();
      ewald.computeEnergy();
      double kEnergyRef = ewald.kSpaceEnergy()/nAtom;
      pair.unsetEnergy();
      pair.computeEnergy();
      double rEnergyRef = spme.rSpaceEnergy()/nAtom;
      double energyRef = kEnergyRef + rEnergyRef;

      // Compute reference force, well converged system
      sim.system().setZeroForces();
      ewald.addForces();
      pair.addForces();
      storeForces(forces_);

      // Loop over values of alpha
      double dAlpha = (alphaMax - alphaMin)/double(n);
      double alpha, kEnergy, rEnergy, energy, fError;
      double rPressure, kPressure, pressure, dPressure;
      for (int i = 0; i <= n; ++i) {
         alpha = alphaMin + dAlpha*i;
         spme.set("alpha", alpha);

         // Energy
         spme.unsetEnergy();
         kEnergy = spme.kSpaceEnergy()/nAtom;
         pair.unsetEnergy();
         rEnergy = spme.rSpaceEnergy()/nAtom;
         energy = kEnergy + rEnergy;

         spme.unsetStress();
         spme.computeStress();
         kPressure = spme.kSpaceStress().trace()*volume/(3.0*nAtom);
         pair.unsetStress();
         pair.computeStress();
         rPressure = spme.rSpaceStress().trace()*volume/(3.0*nAtom);
         pressure = kPressure + rPressure;
         dPressure = pressure - energyRef/3.0;

         sim.system().setZeroForces();
         spme.addForces();
         pair.addForces();
         fError = computeForceError(forces_);

         std::cout << Dbl(spme.get("alpha"), 10)
                   << "  " << Dbl(energy - energyRef, 15) 
                   << "  " << Dbl(fError, 15) 
                   << "  " << Dbl(dPressure, 15) 
                   << std::endl;
      }

      // Reset to default
      spme.set("alpha", alpha_);
   }
 
   MdSimulation& simulation()
   {  return sim; }

private:

   MdSimulation sim;
   MdEwaldPotential ewald;

   DArray<int> capacities_; 
   DArray<double> diameters_; 
   DArray<Vector> forces_;
   double alpha_;
   int nAtom_;

};

int main(int argc, char* argv[]) 
{
   SpmeTest test;

   test.readParam("in/param.spme");
   test.generateConfig();
   test.compareKSpace();
   test.varyAlpha(0.5, 1.2, 14);

}
