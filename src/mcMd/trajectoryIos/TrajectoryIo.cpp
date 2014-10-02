#ifndef MCMD_TRAJECTORY_IO_CPP
#define MCMD_TRAJECTORY_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryIo.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor. 
   */
   TrajectoryIo::TrajectoryIo(System& system)
    : boundaryPtr_(&system.boundary()),
      systemPtr_(&system),
      simulationPtr_(&system.simulation()),
      atomCapacity_(0)
   {}

   /* 
   * Destructor.   
   */
   TrajectoryIo::~TrajectoryIo() 
   {}

   /*
   * Read the trajectory file header.
   */
   void TrajectoryIo::readHeader(std::fstream &file)
   {
      int nSpecies = simulation().nSpecies();
      int speciesCapacity = 0;
      int iSpec,iMol;
      Species* speciesPtr;
      Molecule* molPtr;

      // Add molecules to system
      atomCapacity_ = 0;
      for (iSpec = 0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         speciesCapacity = speciesPtr->capacity();

         for (iMol = 0; iMol < speciesCapacity; ++iMol) {
            molPtr = &(speciesPtr->reservoir().pop());
            system().addMolecule(*molPtr);
         }
         atomCapacity_ += speciesCapacity*speciesPtr->nAtom();
      }

   }

   /*
   * Read a single frame.
   */
   bool TrajectoryIo::readFrame(std::fstream& file)
   {
       UTIL_THROW("This TrajectoryIo class does not implement a readFrame() method.");
       return false;
   }

   /*
   * Write the trajectory file header.
   */
   void TrajectoryIo::writeHeader(std::fstream &file)
   {
      UTIL_THROW("This TrajectoryIo class does not implement a writeHeader() method.");
   }

   /*
   * Write a single frame.
   */
   void TrajectoryIo::writeFrame(std::fstream& file)
   {
       UTIL_THROW("This TrajectoryIo class does not implement a writeFrame() method.");
   }

} 
#endif
