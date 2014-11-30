/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"

namespace McMd
{

   using namespace Util;

   class Species;
   class System;

   // Inline member functions

   /// Constructor.
   Molecule::Molecule()
    : speciesPtr_(0),
      systemPtr_(0),
      firstAtomPtr_(0),
      #ifdef INTER_BOND
      firstBondPtr_(0),
      #endif
      #ifdef INTER_ANGLE
      firstAnglePtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      firstDihedralPtr_(0),
      #endif
      nAtom_(0),
      #ifdef INTER_BOND
      nBond_(0),
      #endif
      #ifdef INTER_ANGLE
      nAngle_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      nDihedral_(0),
      #endif
      id_(NullIndex)
   {}

   #if 0
   /*
   * Set global index.
   */
   void Molecule::setId(int id)
   {  id_ = id; }
   #endif

   /*
   * Mark all atoms and groups in this molecule as active.
   */
   void Molecule::setIsActive()
   {
      int i;
      for (i = 0; i < nAtom_; ++i) {
         (firstAtomPtr_ + i)->setIsActive(true);
      }
      #ifdef INTER_BOND
      for (i = 0; i < nBond_; ++i) {
         (firstBondPtr_ + i)->activate();
      }
      #endif
      #ifdef INTER_ANGLE
      for (i = 0; i < nAngle_; ++i) {
         (firstAnglePtr_ + i)->activate();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      for (i = 0; i < nDihedral_; ++i) {
         (firstDihedralPtr_ + i)->activate();
      }
      #endif
   }
} 
