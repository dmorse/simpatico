/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   TrajectoryReader::TrajectoryReader(System& system)
    : nAtomTotal_(0),
      boundaryPtr_(&system.boundary()),
      systemPtr_(&system),
      simulationPtr_(&system.simulation())
   {}

   /* 
   * Destructor.   
   */
   TrajectoryReader::~TrajectoryReader() 
   {}

   /*
   * Add all molecules and set nAtomTotal_.
   */
   void TrajectoryReader::addMolecules()
   {
      int nSpecies = simulation().nSpecies();
      int speciesCapacity = 0;
      int iSpec,iMol;
      Species* speciesPtr;
      Molecule* molPtr;

      // Add molecules to system
      nAtomTotal_ = 0;
      for (iSpec = 0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         speciesCapacity = speciesPtr->capacity();

         for (iMol = 0; iMol < speciesCapacity; ++iMol) {
            molPtr = &(simulation().getMolecule(iSpec));
            system().addMolecule(*molPtr);
         }
         nAtomTotal_ += speciesCapacity*speciesPtr->nAtom();
      }

   }

} 
