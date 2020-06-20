/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpTrajectoryWriter.h"
#include <simp/species/Species.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/space/Vector.h>

#include <climits>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SmpTrajectoryWriter::SmpTrajectoryWriter(MdSystem& system)
    : TrajectoryWriter(system, true),
      nAtom_(0)
   {  setClassName("SmpTrajectoryWriter"); }

   /*
   * Destructor.
   */
   SmpTrajectoryWriter::~SmpTrajectoryWriter()
   {}

   /*
   * Write data that should appear once, at beginning of the file.
   */
   void SmpTrajectoryWriter::writeHeader()
   {
      BinaryFileOArchive ar(outputFile_);

      // Write Atom Types
      bool hasAtomTypes = true;
      ar << hasAtomTypes;
      bool hasMass = true;
      bool hasCharge = false;
      ar << hasMass;
      ar << hasCharge;
      int nAtomType = simulation().nAtomType();
      ar << nAtomType;
      double mass = 1.0;
      for (int i = 0; i < nAtomType; ++i) {
         ar << mass;
      }

      // Write species information
      bool hasSpecies = true;
      ar << hasSpecies;
      int nSpecies = simulation().nSpecies();
      int nMolecule;
      ar << nSpecies;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         ar << system().nMolecule(iSpecies);
         ar << simulation().species(i);
      }

      // Count and output nAtom_, total number of atoms
      nAtom_ = 0;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         nMolecule = system().nMolecule(iSpecies);
         nAtom_ += nMolecule*simulation().species(iSpecies).nAtom();
      }
      ar << nAtom_;
   }

   void SmpTrajectoryWriter::writeFrame(long iStep)
   {

      BinaryFileOArchive ar(outputFile_);

      ar << iStep;
      ar << boundary();

      // Define and write atom format
      bool isOrdered = true;
      bool hasAtomId = false;
      bool hasAtomContext = false;
      bool hasAtomTypeId = false;
      bool hasAtomVelocity = false;
      bool hasAtomShift = false;
      ar << isOrdered;
      ar << hasAtomId;
      ar << hasAtomContext;
      ar << hasAtomTypeId;
      ar << hasAtomVelocity;
      ar << hasAtomShift;

      // Write total number of atoms
      ar << nAtom_;

      // Write all atoms (ordered)
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      Vector r, v;
      int atomId, j;
      int iAtom, iMolecule, typeId;
      int nSpecies = simulation().nSpecies();
      unsigned int ir;
      atomId = 0; // Global atom index
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         iMolecule = 0;
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            iAtom = 0; // Intramolecule atom index
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               if (hasAtomId) {
                  ar << atomId;
               }
               if (hasAtomContext) {
                  ar << iSpecies;
                  ar << iMolecule;
                  ar << iAtom;
               }
               if (hasAtomTypeId) {
                  typeId = atomIter->typeId();
                  ar << typeId;
               }
               boundary().transformCartToGen(atomIter->position(), r);
               for (j = 0; j < Dimension; ++j) {
                  if (r[j] >= 1.0) r[j] -= 1.0;
                  if (r[j] <  0.0) r[j] += 1.0;
                  ir = floor( UINT_MAX*r[j] + r[j] + 0.5 );
                  ar << ir;
               }
               if (hasAtomVelocity) {
                  v = atomIter->velocity();
                  ar << v;
               }
               // Shift
               ++iAtom;
               ++atomId;
            } // atom loop
            ++iMolecule;
         } // molecule loop
      } // species loop

   }

}
