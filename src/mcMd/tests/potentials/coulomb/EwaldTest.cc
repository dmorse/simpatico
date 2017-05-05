
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

   MdSimulation sim;

public:


   void setup() 
   {
      // Read parameter file
      std::ifstream in;
      sim.fileMaster().openInputFile("in/param", in);
      sim.readParam(in);
      in >> sim.system().boundary();
      in.close();
      std::cout << "Finished reading param file" << std::endl;

      // Generate an initial configuration
      DArray<int> capacities; 
      int nSpecies = sim.nSpecies();
      capacities.allocate(nSpecies);
      DArray<double> diameters; 
      int nAtomType = sim.nAtomType();
      diameters.allocate(nAtomType);
      for (int i = 0; i < nSpecies; ++i) {
         capacities[i] = sim.species(i).capacity();
      }
      for (int i = 0; i < nAtomType; ++i) {
         diameters[i] = 0.2;
      }
      std::cout << "Begin generating molecules" << std::endl;
      sim.system().generateMolecules(capacities, diameters);
      std::cout << "Finished generating molecules" << std::endl;

      double alphaMin = 0.5;
      double dAlpha   = 0.05;
      double alpha, kEnergy, rEnergy, energy;
      MdCoulombPotential& coulomb = sim.system().coulombPotential();
      MdPairPotential& pair = sim.system().pairPotential();
      for (int i = 0; i < 20; ++i) {
         alpha = alphaMin + dAlpha*i;
         coulomb.set("alpha", alpha);
         coulomb.unsetEnergy();
         pair.unsetEnergy();
         kEnergy = coulomb.kSpaceEnergy();
         rEnergy = coulomb.rSpaceEnergy();
         energy = kEnergy + rEnergy;
         std::cout << Dbl(coulomb.get("alpha"),15)
                   << "  " << Dbl(kEnergy, 15)
                   << "  " << Dbl(rEnergy, 15)
                   << "  " << energy
                   << std::endl;
      }
   }

};

int main(int argc, char* argv[]) 
{

   EwaldTest test;

   test.setup();
   

}
