#ifdef SIMP_EXTERNAL
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalEnergyAverage.h"        // class header

#include <util/misc/FileMaster.h>  
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/external/ExternalPotential.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McExternalEnergyAverage::McExternalEnergyAverage(McSystem& system)
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McExternalEnergyAverage"); }
 
   /*
   * Evaluate external energy, and add to ensemble. 
   */
   void McExternalEnergyAverage::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;

      double energy = 0.0;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
             for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                 energy += system().externalPotential().energy(atomIter->position(), atomIter->typeId());
             }
         }
      }
      accumulator_.sample(energy, outputFile_);
   }

}
#endif 
