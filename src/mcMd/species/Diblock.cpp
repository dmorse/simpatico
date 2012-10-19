#ifndef MCMD_DIBLOCK_CPP
#define MCMD_DIBLOCK_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
    #ifdef INTER_ANGLE
    , angleType_(NullIndex)
    #endif
    #ifdef INTER_DIHEDRAL
    , dihedralType_(NullIndex)
    #endif
   { 
      setClassName("Diblock"); 
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

      #if INTER_ANGLE
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

      #if INTER_DIHEDRALS
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
   * Load from Serializable::IArchive.
   */
   void Diblock::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadCArray<int>(ar, "blockLengths", blockLengths_, 2);
      for (int i=0; i < 2; ++i) {
         if (blockLengths_[i] < 1) {
            UTIL_THROW("Invalid blockLength for diblock.");
         }
      }
      nAtom_ = blockLengths_[0] + blockLengths_[1];
      loadCArray<int>(ar, "atomTypes", atomTypes_, 2);

      nBond_  = nAtom_ - 1;
      loadParameter<int>(ar,"bondType", bondType_);
      #ifdef INTER_ANGLE
      loadParameter<int>(ar,"hasAngles", hasAngles_);
      if (hasAngles_) {
         nAngle_ = nBond_ - 1;
         if (nAngle_ > 0) {
            loadParameter<int>(ar,"angleType", angleType_);
         }
      } else {
         nAngle_ = 0;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      loadParameter<int>(ar,"hasDihedrals", hasDihedrals_);
      if (hasDihedrals_) {
         if (nAtom_ > 3) {
            nDihedral_ = nAtom_ - 3;
         } else {
            nDihedral_ = 0;
         }
         if (nDihedral_ > 0) {
            loadParameter<int>(ar, "dihedralType", dihedralType_);
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
   * Save internal state to an archive.
   */
   void Diblock::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar.pack(blockLengths_, 2);
      ar.pack(atomTypes_, 2);
      ar << bondType_;
      #ifdef INTER_ANGLE
      ar << hasAngles_;
      if (hasAngles_ && nAngle_ > 0) {
         ar << angleType_;
      } 
      #endif
      #ifdef INTER_DIHEDRAL
      ar << hasDihedrals_;
      if (hasDihedrals_ && nDihedral_ > 0) {
         ar << dihedralType_;
      } 
      #endif
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

   #ifdef INTER_ANGLE
   /* 
   * Return angleType_ for every angle.
   */
   int Diblock::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /* 
   * Return dihedralType_ for every dihedral.
   */
   int Diblock::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif

} 
#endif
