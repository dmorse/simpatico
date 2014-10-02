#ifndef MCMD_LAMMPS_DUMP_IO_CPP
#define MCMD_LAMMPS_DUMP_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpIo.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpIo::LammpsDumpIo(System &system)
   : TrajectoryIo(system)
   {}

   /*
   * Destructor.
   */
   LammpsDumpIo::~LammpsDumpIo()
   {}

   void LammpsDumpIo::readHeader(std::fstream &file)
   {}

   void LammpsDumpIo::readFrame(std::fstream &file)
   {

      bool notEnd;
      std::stringstream line;

      // Read ITEM: TIMESTEP
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: TIMESTEP");
      }
      checkString(line, "ITEM:");
      checkString(line, "TIMESTEP");
      int step;
      file >> step;
 
      // Read ITEM: NUMBER OF ATOMS
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: NUMBER OF ATOMS");
      }
      checkString(line, "ITEM:");
      checkString(line, "NUMBER");
      checkString(line, "OF");
      checkString(line, "ATOMS");
      int nAtom;
      file >> nAtom;

      if (positions_.isAllocated()) {
         if (nAtom != positions_.capacity()) {
            UTIL_THROW("Inconsistent values of nAtom");
         }
      } else {
         positions_.allocate(nAtom);
      }

      // Read ITEM: BOX
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: BOX");
      }
      checkString(line, "ITEM:");
      checkString(line, "BOX");
      // Ignore rest of ITEM: BOX line
      Vector min, max, lengths;
      for (int i = 0; i < Dimension; ++i) {
         file >> min[i] >> max[i];
         lengths[i] = max[i] - min[i];
      }

      // Read ITEM: ATOMS 
      notEnd = getNextLine(file, line);
      checkString(line, "ITEM:");
      checkString(line, "ATOMS");
      // Ignore the rest of ITEM: ATOMS  line, for now

      // Load positions into positions_ vector, indexed by id.
      IntVector shift;
      int id, typeId, molId, i, j;
      for (i = 0; i < nAtom; ++i) {

         file >> id;
         id = id - 1; // Convert convention, Lammps -> Simpatico
         file >> typeId;
         file >> molId;
 
         // Read position
         for (j = 0; j < Util::Dimension; ++j) {
            file >> positions_[id][j];
         }
         file >> shift;

         // shift into simulation cell
         boundary().shift(positions_[id]);

      }

      // Load atom positions, assuming ids are ordered by molecule and species
      int iSpecies,iMol;
      Species *speciesPtr;
      Molecule::AtomIterator atomIter;
      Molecule *molPtr;
      id = 0;
      for (iSpecies = 0; iSpecies < simulation().nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         for (iMol = 0; iMol < speciesPtr->capacity(); ++iMol) {
            molPtr = &system().molecule(iSpecies, iMol);
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->position() = positions_[id];
               id++;
            }
         }
      }

      //return true;
   }
}
#endif
