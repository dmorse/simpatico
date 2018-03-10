/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeformCommand.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
//#include <util/format/Str.h>
//#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   DeformCommand::DeformCommand(System& system) 
    : Command("DEFORM_CELL"),
      SystemInterface(system)
   {  setClassName("DeformCommand"); }

   /*
   * Default destructor.
   */
   DeformCommand::~DeformCommand()
   {}

   /*
   * Execute command to deform unit cell.
   */
   void DeformCommand::execute(std::istream& in)
   {
 
      // Transform atomic positions to generalized [0,1] coordinates
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

      // Transform positions back to Cartesian
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

      // Rebuild cell and/or pair list.
      reneighbor();

   }

}
