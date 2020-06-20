/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpTrajectoryReader.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/species/SpeciesFinder.h>
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
   SmpTrajectoryReader::SmpTrajectoryReader(System &system)
   : TrajectoryReader(system)
   {}

   /*
   * Destructor.
   */
   SmpTrajectoryReader::~SmpTrajectoryReader()
   {}

   /*
   * Open trajectory file and setup to read.
   */
   void SmpTrajectoryReader::open(std::string filename)
   {
      // Open trajectory file for reading
      simulation().fileMaster().openInputFile(filename, file_);

      // Assocate a binary file input archive with file_
      BinaryFileIArchive ar(file_);

      // Read atom type data (optional)
      bool hasAtomTypes;
      ar >> hasAtomTypes;
      if (hasAtomTypes) {
         double mass, charge;
         int nAtomType = simulation().nAtomType();
         bool hasMass, hasCharge;
         ar >> hasMass;
         ar >> hasCharge;
         ar >> nAtomType;
         for (int i = 0; i < nAtomType; ++i) {
             if (hasMass) ar >> mass;
             if (hasCharge) ar >> charge;
         }
      }

      // Read and check molecular species data (required)
      bool hasSpecies;
      ar >> hasSpecies;
      UTIL_CHECK(hasSpecies);
      int nSpecies;
      ar >> nSpecies;
      UTIL_CHECK(nSpecies = simulation().nSpecies());
      nMolecules_.allocate(nSpecies);
      Species  speciesIn;  // species read from file
      Species* speciesPtr; // species in parent simulation
      for (int i = 0; i < nSpecies; ++i) {
         ar >> nMolecules_[i];
         ar >> speciesIn;
         speciesPtr = &simulation().species(i);
         UTIL_CHECK(speciesIn.nAtom() == speciesPtr->nAtom());
         UTIL_CHECK(nMolecules_[i] == speciesPtr->capacity());
         UTIL_CHECK(nMolecules_[i] <= speciesIn.capacity());
      }

      // Add all molecules to system, set nAtomTotal_
      addMolecules();

      // Read and check total number of atoms
      int nAtom;
      ar >> nAtom;
      UTIL_CHECK(nAtom == nAtomTotal_);
     
   }

   /*
   * Read frame, return false if end-of-file
   */
   bool SmpTrajectoryReader::readFrame()
   {
      // Associate binary file input archive with file_
      BinaryFileIArchive ar(file_);

      // Attempt to read iStep
      long iStep = -1;  
      ar >> iStep;

      // Return false if read failed, indicating end of file.
      if (file_.eof()) {
         return false;
      }

      // Read boundary dimensions (period unit cell)
      ar >> boundary();

      // Read atom data format
      bool isOrdered;
      bool hasAtomId;
      bool hasAtomContext;
      bool hasAtomTypeId;
      bool hasAtomVelocity;
      bool hasAtomShift;
      ar >> isOrdered;        // Are atom ids in consecutive order
      ar >> hasAtomId;        // Is the global atom id included
      ar >> hasAtomContext;   // Is molecular context included
      ar >> hasAtomTypeId;    // Is the atom type Id included
      ar >> hasAtomVelocity;  // Is the velocity included
      ar >> hasAtomShift;     // Is an atom shift included
      if (!isOrdered) UTIL_CHECK(hasAtomId);

      // Set up SpeciesFinder 
      Species *speciesPtr;
      SpeciesFinder finder;
      SpeciesFinder::Part context;
      int nSpecies = simulation().nSpecies();
      finder.allocate(nSpecies);
      int nMolecule, nAtomMol;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         nMolecule = nMolecules_[iSpecies];
         nAtomMol  = speciesPtr->nAtom();
         finder.setSpecies(iSpecies, nMolecule, nAtomMol);
      }
      finder.initialize();

      // Read total number of atoms
      int nAtom;
      ar >> nAtom;
      UTIL_CHECK(nAtom == nAtomTotal_);

      // Read all atoms 
      Molecule *molPtr;
      Atom *atomPtr;
      Vector r, v;
      double h = 1.0/(double(UINT_MAX) + 1.0);
      int atomId, iSpecies, iAtom, iMolecule, typeId, j;
      unsigned int ir;
      for (int i = 0; i < nAtomTotal_; ++i) {
         if (hasAtomId) {
            ar >> atomId;
            if (isOrdered) UTIL_CHECK(atomId == i);
         } else {
            atomId = i;
         }
         finder.findPart(atomId, context);
         if (hasAtomContext) {
            ar >> iSpecies;
            UTIL_CHECK(iSpecies == context.speciesId);
            ar >> iMolecule;
            UTIL_CHECK(iMolecule == context.moleculeId);
            ar >> iAtom;
            UTIL_CHECK(iAtom == context.partId);
         } else {
            iSpecies = context.speciesId;
            iMolecule = context.moleculeId;
            iAtom = context.partId;
         }
         molPtr = &(system().molecule(iSpecies, iMolecule));
         atomPtr = &(molPtr->atom(iAtom));
         if (hasAtomTypeId) {
            ar >> typeId;
            atomPtr->setTypeId(typeId);
         }
         // Read position, scaled componenets in uint form
         for (j = 0; j < Dimension; ++j) {
            ar >> ir;
            r[j] = ir*h;
         }
         boundary().transformGenToCart(r, atomPtr->position());
         if (hasAtomVelocity) {
            ar >> atomPtr->velocity();
         }
         // Read shift, if any
      }

      return true;
   }

   /*
   * Close trajectory file.
   */
   void SmpTrajectoryReader::close()
   {  file_.close(); }

}
