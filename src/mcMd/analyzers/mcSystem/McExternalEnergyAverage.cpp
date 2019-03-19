#ifdef SIMP_EXTERNAL
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalEnergyAverage.h"     
#include <mcMd/mcSimulation/McSystem.h>
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
    : McAverageAnalyzer(system)
   {  setClassName("McExternalEnergyAverage"); }

   /*
   * Evaluate external energy
   */
   void McExternalEnergyAverage::compute() 
   {
      double value_ = 0.0;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      //ExternalPotential& potential = system().externalPotential();
      for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
             for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                 value_ += system().externalPotential().energy(atomIter->position(), atomIter->typeId());
                 //value_ += potential().energy(atomIter->position(), atomIter->typeId());
             }
         }
      }
   }

}
#endif 
