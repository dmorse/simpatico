
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/species/Species.h>

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
         diameters[i] = 1.0;
      }
      std::cout << "Begin generating molecules" << std::endl;
      sim.system().generateMolecules(capacities, diameters);
      std::cout << "Finished generating molecules" << std::endl;
   }

};

int main(int argc, char* argv[]) 
{

   EwaldTest test;

   test.setup();

}
