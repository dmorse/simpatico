#ifndef DIBLOCK_CPP
#define DIBLOCK_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diblock.h"
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   Diblock::Diblock()
    : Linear()
    #ifdef MCMD_ANGLE
    , angleType_(NullIndex)
    #endif
    #ifdef MCMD_DIHEDRAL
    , dihedralType_(NullIndex)
    #endif
   { 
      for (int i=0; i < 2; ++i) {
         atomTypes_[i] = NullIndex;  
         blockLengths_[i] = 0; 
      }
   }
   
   /* 
   * Read blockLengths_ and atomTypes_ for each block
   */
   void Diblock::readSpeciesParam(std::istream &in)
   {
      readCArray<int>(in,"blockLengths", blockLengths_, 2);
      for (int i=0; i < 2; ++i) {
         if (blockLengths_[i] < 1) {
            UTIL_THROW("Invalid blockLength for diblock.");
         }
      }
      nAtom_ = blockLengths_[0] + blockLengths_[1];
      readCArray<int>(in,"atomTypes", atomTypes_, 2);

      nBond_ = nAtom_ - 1;
      read<int>(in, "bondType", bondType_);

      #if MCMD_ANGLE
      read<int>(in, "hasAngle", hasAngles_);
      if (nAtom_ <= 3) {
         UTIL_THROW("Error: Cannot have dihedrals with nAtom <= 4");
      }
      if (hasAngles_) {
         nAngle_ = nAtom_ - 2;
         read<int>(in, "angleType", angleType_);
      } else {
         nAngle_ = 0;
      }
      #endif

      #if MCMD_DIHEDRALS
      read<int>(in, "hasDihedrals", hasDihedrals_);
      if (nAtom_ <= 4) {
         UTIL_THROW("Error: Cannot have dihedrals with nAtom <= 4");
      }
      if (hasDihedrals_) {
         nDihedral_ = nAtom_ - 3;
         read<int>(in, "dihedralType", dihedralType_);
      } else {
         nDihedral_ = 0;
      }
      #endif

      buildLinear();
   }
 
   /* 
   * Return atomTypes_ for every atom.
   */
   int Diblock::calculateAtomTypeId(int index) const
   { 
      if (index < blockLengths_[0]) {
         return atomTypes_[0];
      } else {
         return atomTypes_[1];
      }
   }
   
   /* 
   * Return bondType_ for every bond.
   */
   int Diblock::calculateBondTypeId(int molBondId) const
   { return bondType_; }

   #ifdef MCMD_ANGLE
   /* 
   * Return angleType_ for every angle.
   */
   int Diblock::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef MCMD_DIHEDRAL
   /* 
   * Return dihedralType_ for every dihedral.
   */
   int Diblock::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif

} 
#endif
