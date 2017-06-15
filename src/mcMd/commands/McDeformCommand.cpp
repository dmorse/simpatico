/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McDeformCommand.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <util/format/Str.h>
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   McDeformCommand::McDeformCommand(McSystem& system) 
    : McSystemInterface(system)
   {}

   /*
   * Default destructor.
   */
   McDeformCommand::~McDeformCommand()
   {}

   /*
   * Read and execute command to deform unit cell.
   */
   bool 
   McDeformCommand::readCommand(std::string const & name, std::istream& in)
   { 
      if (name != "DEFORM_CELL") {
 
         #if 0
         // Read in configuration from file
         std::string filename;
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         std::ifstream inputFile;
         system().fileMaster().openInputFile(filename, inputFile);
         system().readConfig(inputFile);
         inputFile.close();
         #endif

         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;
         Vector cartPosition, genPosition;
         int nSpecies = simulation().nSpecies();
         for (int iSpec=0; iSpec < nSpecies; ++iSpec) {
            begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  cartPosition = atomIter->position();
                  boundary().transformCartToGen(cartPosition, genPosition);
                  atomIter->position() = genPosition;
               }
            }
         }

         // Read new system boundary
         in >> boundary();
         Log::file() << "  " << system().boundary();
         Log::file() << std::endl;

         // Convert positions back to Cartesian
         for (int iSpec=0; iSpec < nSpecies; ++iSpec) {
            begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  genPosition = atomIter->position();
                  boundary().transformGenToCart(genPosition, cartPosition);
                  atomIter->position() = cartPosition;
               }
            }
         }

         #ifndef SIMP_NOPAIR 
         // Generate cell list
         pairPotential().buildCellList();
         #endif

         #if 0
         // Write out configuration to file
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         std::ofstream outputFile;
         system().fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();
         #endif

         // Command name recognized, successful completion.
         return true;

      } else {

         // if command name not recognized
         return false;

      }
   }

}
