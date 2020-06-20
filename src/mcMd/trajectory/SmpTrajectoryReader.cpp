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
      // Open trajectory file
      simulation().fileMaster().openInputFile(filename, file_);
      BinaryFileIArchive ar(file_);

      // Read atom type data
      bool hasAtomTypes;
      bool hasMass;
      bool hasCharge;
      ar >> hasAtomTypes;
      ar >> hasMass;
      ar >> hasCharge;
      int nAtomType = simulation().nAtomType();
      ar >> nAtomType;
      double mass, charge;
      for (int i = 0; i < nAtomType; ++i) {
          if (hasMass) ar >> mass;
          if (hasCharge) ar >> charge;
      }

      #if 0
      bool hasSpecies = true;
      ar >> hasSpecies;
      ar >> nSpecies;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      }
      #endif

      int nAtom;
      ar >> nAtom;
      
      // Add all molecules to system and check consistency of nAtom.
      addMolecules();
      UTIL_CHECK(nAtom == nAtomTotal_);
     
   }

   /*
   * Read frame, return false if end-of-file
   */
   bool SmpTrajectoryReader::readFrame()
   {
      // Preconditions

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

      // Read atom format
      bool isOrdered = true;
      bool hasAtomId = false;
      bool hasAtomContext = false;
      bool hasAtomTypeId = false;
      bool hasAtomVelocity = false;
      bool hasAtomShift = false;
      ar >> isOrdered;
      ar >> hasAtomId;
      ar >> hasAtomContext;
      ar >> hasAtomTypeId;
      ar >> hasAtomVelocity;
      ar >> hasAtomShift;

      // Set up SpeciesFinder 
      Species *speciesPtr;
      SpeciesFinder finder;
      SpeciesFinder::Part context;
      int nSpecies = simulation().nSpecies();
      finder.allocate(nSpecies);
      int nMolecule, nAtomMol;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         nMolecule = speciesPtr->capacity();
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
         }
         molPtr = &(system().molecule(context.speciesId, context.moleculeId));
         atomPtr = &(molPtr->atom(context.partId));
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
