#ifndef HOMOPOLYMER_CPP
#define HOMOPOLYMER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Homopolymer.h"

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   Homopolymer::Homopolymer()
    : Linear(),
      atomType_(NullIndex),
      bondType_(NullIndex)
      #ifdef MCMD_ANGLE
      , angleType_(NullIndex)
      #endif
      #ifdef MCMD_DIHEDRAL
      , dihedralType_(NullIndex)
      #endif
   {} 
   
   /* 
   * Read nAtom and type.
   */
   void Homopolymer::readSpeciesParam(std::istream &in)
   {
      read<int>(in,"nAtom", nAtom_);
      read<int>(in,"atomType", atomType_);
      nBond_  = nAtom_ - 1;
      read<int>(in,"bondType", bondType_);
      #ifdef MCMD_ANGLE
      read<int>(in,"hasAngles", hasAngles_);
      if (hasAngles_) {
         nAngle_ = nBond_ - 1;
         if (nAngle_ > 0) {
            read<int>(in,"angleType", angleType_);
         }
      } else {
         nAngle_ = 0;
      }
      #endif
      #ifdef MCMD_DIHEDRAL
      read<int>(in,"hasDihedrals", hasDihedrals_);
      if (hasDihedrals_) {
         if (nAtom_ > 3) {
            nDihedral_ = nAtom_ - 3;
         } else {
            nDihedral_ = 0;
         }
         if (nDihedral_ > 0) {
            read<int>(in,"dihedralType", dihedralType_);
         } else {
            nDihedral_ = 0;
         }
      } else {
         nDihedral_ = 0;
      }
      #endif

      buildLinear();
   }

   /* 
   * Return atomType_ for every atom.
   */
   int Homopolymer::calculateAtomTypeId(int index) const
   { return atomType_; }

   /* 
   * Return bondType_ for every bond.
   */
   int Homopolymer::calculateBondTypeId(int index) const
   { return bondType_; }

   #ifdef MCMD_ANGLE
   /* 
   * Return angleType_ for every angle.
   */
   int Homopolymer::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef MCMD_DIHEDRAL
   /* 
   * Return dihedralType_ for every dihedral.
   */
   int Homopolymer::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif
 
} 
#endif
