/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpReader.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
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
   LammpsDumpReader::LammpsDumpReader(System &system)
   : TrajectoryReader(system)
   {}

   /*
   * Destructor.
   */
   LammpsDumpReader::~LammpsDumpReader()
   {}

   /*
   * Open file and setup memory.
   */
   void LammpsDumpReader::open(std::string filename)
   {
      // Open trajectory file
      simulation().fileMaster().openInputFile(filename, file_);

      // Set nAtomTotal_ and add all molecules to system
      addMolecules();

      // Allocate private array of atomic positions_
      if (!positions_.isAllocated()) {
         positions_.allocate(nAtomTotal_);
      } else {
         if (nAtomTotal_ != positions_.capacity()) {
            UTIL_THROW("Inconsistent values of atom capacity");
         }
      }
   }

   /*
   * Read frame, return false if end-of-file
   */
   bool LammpsDumpReader::readFrame()
   {
      // Preconditions
      if (!positions_.isAllocated()) {
         UTIL_THROW("positions_ array is not allocated");
      }

      bool notEnd;
      std::stringstream line;

      // Attempt to read first line
      notEnd = getNextLine(file_, line);
      if (!notEnd) {
         return false;
      }

      // Process ITEM: TIMESTEP
      checkString(line, "ITEM:");
      checkString(line, "TIMESTEP");
      int step;
      file_ >> step;
 
      // Read ITEM: NUMBER OF ATOMS
      notEnd = getNextLine(file_, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: NUMBER OF ATOMS");
      }
      checkString(line, "ITEM:");
      checkString(line, "NUMBER");
      checkString(line, "OF");
      checkString(line, "ATOMS");
      int nAtom;
      file_ >> nAtom;
      if (nAtom != nAtomTotal_) {
         UTIL_THROW("Inconsistent values: nAtom != nAtomTotal_");
      }

      // Read ITEM: BOX
      notEnd = getNextLine(file_, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: BOX");
      }
      checkString(line, "ITEM:");
      checkString(line, "BOX");
      // Ignore rest of ITEM: BOX line
      Vector min, max, lengths;
      for (int i = 0; i < Dimension; ++i) {
         file_ >> min[i] >> max[i];
         lengths[i] = max[i] - min[i];
      }
      boundary().setOrthorhombic(lengths);

      // Read ITEM: ATOMS 
      notEnd = getNextLine(file_, line);
      checkString(line, "ITEM:");
      checkString(line, "ATOMS");
      // Ignore the rest of ITEM: ATOMS  line, for now

      // Load positions into positions_ vector, indexed by id.
      IntVector shift;
      int id, typeId, molId, i, j;
      for (i = 0; i < nAtom; ++i) {

         file_ >> id;
         id = id - 1; // Convert convention, Lammps -> Simpatico
         file_ >> typeId;
         file_ >> molId;
 
         // Read position
         for (j = 0; j < Util::Dimension; ++j) {
            file_ >> positions_[id][j];
         }
         file_ >> shift;

         // shift into simulation cell
         boundary().shift(positions_[id]);

      }

      // Assign atom positions, assuming ordered atom ids 
      int iSpecies, iMol;
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

      return true;
   }

   /*
   * Close trajectory file.
   */
   void LammpsDumpReader::close()
   {  file_.close();}

}
