/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryWriter.h"
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
   DdMdTrajectoryWriter::DdMdTrajectoryWriter(MdSystem& system)
    : TrajectoryWriter(system, true),
      nAtomTot_(0)
   {  setClassName("DdMdTrajectoryWriter"); }

   /*
   * Destructor.
   */
   DdMdTrajectoryWriter::~DdMdTrajectoryWriter()
   {}

   /*
   * Write data that should appear once, at beginning of the file.
   */
   void DdMdTrajectoryWriter::writeHeader()
   {
      // Count and output total number of atoms, nAtomTot_
      int nMolecule;
      int nSpecies = simulation().nSpecies();
      nAtomTot_ = 0;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         nMolecule = system().nMolecule(iSpecies);
         nAtomTot_ += nMolecule*simulation().species(iSpecies).nAtom();
      }

      BinaryFileOArchive ar(outputFile_);
      ar << nAtomTot_;
   }

   void DdMdTrajectoryWriter::writeFrame(long iStep)
   {

      BinaryFileOArchive ar(outputFile_);

      ar << iStep;
      ar << boundary();

      // Write ATOMS block
      // ar << endl << "ordered";
      // ar << endl << "format imtpv";
      // ar << endl << "nAtom " << nAtomTot;

      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      Vector r;
      int atomId, j;
      // int iAtom, iMolecule, typeId;
      int nSpecies = simulation().nSpecies();
      unsigned int ir;
      atomId = 0; // Global atom index
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         // iMolecule = 0;
         system().begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            // iAtom = 0; // Intramolecule atom index
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               ar << atomId;
               // ar << iSpecies;
               // ar << iMolecule;
               // ar << iAtom;
               // typeId = atomIter->typeId();
               // ar << typeId;
               boundary().transformCartToGen(atomIter->position(), r);
               for (j = 0; j < Dimension; ++j) {
                  if (r[j] >= 1.0) r[j] -= 1.0;
                  if (r[j] <  0.0) r[j] += 1.0;
                  ir = floor( UINT_MAX*r[j] + r[j] + 0.5 );
                  ar << ir;
               }
               // ar << atomIter->velocity();
               // ++iAtom;
               ++atomId;
            } // atom loop
            // ++iMolecule;
         } // molecule loop
      } // species loop

   }

}
