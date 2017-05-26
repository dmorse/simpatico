/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryReader.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/misc/ioUtil.h>

#include <climits>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdTrajectoryReader::DdMdTrajectoryReader(System &system)
   : TrajectoryReader(system)
   {}

   /*
   * Destructor.
   */
   DdMdTrajectoryReader::~DdMdTrajectoryReader()
   {}

   /*
   * Open trajectory file and setup to read.
   */
   void DdMdTrajectoryReader::open(std::string filename)
   {
      // Open trajectory file
      simulation().fileMaster().openInputFile(filename, file_);

      BinaryFileIArchive ar(file_);
      int nAtom;
      ar >> nAtom;
      //std::cout << nAtom << std::endl;
      
      // Add all molecules to system and check consistency of nAtom.
      addMolecules();
      if (nAtom != nAtomTotal_) {
         UTIL_THROW("Inconsistent values: nAtom != nAtomTotal_");
      }
     
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
   bool DdMdTrajectoryReader::readFrame()
   {
      // Preconditions
      if (!positions_.isAllocated()) {
         UTIL_THROW("positions_ array is not allocated");
      }

      BinaryFileIArchive ar(file_);

      // Attempt to read iStep
      long iStep = -1;  
      ar >> iStep;

      // Return false if read failed, indicating end of file.
      if (file_.eof()) {
         return false;
      }

      // Read boundary dimensions
      ar >> boundary();

      // Loop over atoms, read atomic positions
      Vector r;
      double h = 1.0/(double(UINT_MAX) + 1.0);
      int id, i, j;
      unsigned int ir;
      for (i = 0; i < nAtomTotal_; ++i) {
         ar >> id;
         for (j = 0; j < Dimension; ++j) {
            ar >> ir;
            r[j] = ir*h;
         }
         boundary().transformGenToCart(r, positions_[id]);
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
   void DdMdTrajectoryReader::close()
   {  file_.close(); }

}
