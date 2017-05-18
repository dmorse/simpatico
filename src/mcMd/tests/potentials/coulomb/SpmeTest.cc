#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
#include <mcMd/potentials/coulomb/MdSpmePotential.h>
#include <mcMd/potentials/coulomb/MdEwaldPotential.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/species/Species.h>
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
      //rSpaceCutoff_ = spme.get("rSpaceCutoff");
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

      sim.system().setZeroForces();
      ewald.addForces();
      storeForces(forces_);

      spme.unsetEnergy();
      spme.makeWaves();
      spme.computeEnergy();
      double kEnergy = spme.kSpaceEnergy();

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


      // Compute properties of a well converged system
      spme.set("alpha", ewald.get("alpha"));
      std::cout << "Ewald reference alpha   = " << ewald.get("alpha") << std::endl;
      std::cout << "Ewald reference epsilon = " << ewald.get("epsilon") << std::endl;
      std::cout << "SPME  reference alpha   = " << spme.get("alpha")  << std::endl;
      std::cout << "SPME  reference epsilon = " << spme.get("epsilon") << std::endl;
      ewald.unsetEnergy();
      ewald.computeEnergy();
      double kEnergyRef = ewald.kSpaceEnergy();
      pair.unsetEnergy();
      pair.computeEnergy();
      double rEnergyRef = spme.rSpaceEnergy();
      double energyRef = kEnergyRef + rEnergyRef;

      //sim.system().setZeroForces();
      //ewald.addForces();
      //storeForces(forces_);

      // Loop over values of alpha
      double dAlpha = (alphaMax - alphaMin)/double(n);
      double alpha, kEnergy, rEnergy, energy;
      // double fError;
      for (int i = 0; i <= n; ++i) {
         alpha = alphaMin + dAlpha*i;
         spme.set("alpha", alpha);

         spme.unsetEnergy();
         spme.computeEnergy();
         kEnergy = spme.kSpaceEnergy();

         pair.unsetEnergy();
         pair.computeEnergy();
         rEnergy = spme.rSpaceEnergy();
         energy = kEnergy + rEnergy;

         // sim.system().setZeroForces();
         // coulomb.addForces();
         // pair.addForces();
         // fError = computeForceError(forces_);

         std::cout << Dbl(spme.get("alpha"), 10)
                   << "  " << Dbl(energy - energyRef, 15) 
                   //<< "  " << Dbl(fError, 15) 
                   << std::endl;
      }

      //Reset to default
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
   double rSpaceCutoff_;
   int nAtom_;

};

int main(int argc, char* argv[]) 
{
   SpmeTest test;

   test.readParam("in/param.spme");
   test.generateConfig();
   test.compareKSpace();
   test.varyAlpha(0.5, 1.0, 5);

}
