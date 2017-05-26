/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Multiblock.h"
#include <util/global.h>

namespace Simp
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   Multiblock::Multiblock()
    : Linear(),
      nBlock_(2)
    #ifdef SIMP_ANGLE
    , angleType_(NullIndex)
    #endif
    #ifdef SIMP_DIHEDRAL
    , dihedralType_(NullIndex)
    #endif
   {  setClassName("Multiblock");
      
   }
   
   /* 
   * Read blockLengths_ and atomTypes_ for each block
   */
   void Multiblock::readSpeciesParam(std::istream &in)
   {
      read<int>(in,"nBlock", nBlock_);
      blockLengths_.allocate(nBlock_);
      readDArray<int>(in,"blockLengths", blockLengths_, nBlock_);
      for (int i=0; i < nBlock_; ++i) {
         if (blockLengths_[i] < 1) {
            UTIL_THROW("Invalid blockLength for Multiblock.");
         }
      }
      atomTypes_.allocate(nBlock_);
      readDArray<int>(in,"atomTypes", atomTypes_, nBlock_);
      for (int i=0; i < nBlock_; ++i) {
         if (atomTypes_[i] < 0) {
            UTIL_THROW("Invalid atomType for Multiblock.");
         }
      }
      blockBegin_.allocate(nBlock_);
      for (int i=0; i < nBlock_; ++i) {
         if (i == 0) {
            blockBegin_[i] = 0;
         } else {
            blockBegin_[i] = blockBegin_[i-1] + blockLengths_[i-1];
         }
      }
      nAtom_ = 0; 
      for (int i=0; i < nBlock_; ++i) {
         nAtom_ = nAtom_ + blockLengths_[i];
      }

      nBond_ = nAtom_ - 1;
      read<int>(in, "bondType", bondType_);

      #ifdef SIMP_ANGLE
      hasAngles_ = 0; // default value
      readOptional<int>(in, "hasAngles", hasAngles_);
      if (hasAngles_) {
         if (nAtom_ < 3) {
            UTIL_THROW("Error: Cannot have angles with nAtom < 3");
         }
         nAngle_ = nAtom_ - 2;
         read<int>(in, "angleType", angleType_);
      } else {
         nAngle_ = 0;
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      hasDihedrals_ = 0; // default value
      readOptional<int>(in, "hasDihedrals", hasDihedrals_);
      if (hasDihedrals_) {
         if (nAtom_ < 4) {
            UTIL_THROW("Error: Cannot have dihedrals with nAtom < 4");
         }
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
   void Multiblock::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar,"nBlock", nBlock_);
      blockLengths_.allocate(nBlock_);
      loadDArray<int>(ar, "blockLengths", blockLengths_, nBlock_);
      for (int i=0; i < nBlock_; ++i) {
         if (blockLengths_[i] < 1) {
            UTIL_THROW("Invalid blockLength for Multiblock.");
         }
      }
      atomTypes_.allocate(nBlock_);
      loadDArray<int>(ar, "atomTypes", atomTypes_, nBlock_);
      for (int i=0; i < nBlock_; ++i) {
         if (atomTypes_[i] < 0) {
            UTIL_THROW("Invalid atomType for Multiblock.");
         }
      }
      blockBegin_.allocate(nBlock_);
      for (int i = 0; i < nBlock_; ++i) {
         if (i == 0) {
            blockBegin_[i] = 0;
         } else {
            blockBegin_[i] = blockBegin_[i-1] + blockLengths_[i-1];
         }
      }
      nAtom_ = 0; 
      for (int i=0; i < nBlock_; ++i) {
         nAtom_ = nAtom_ + blockLengths_[i];
      }

      nBond_ = nAtom_ - 1;
      loadParameter<int>(ar, "bondType", bondType_);

      #ifdef SIMP_ANGLE
      hasAngles_ = 0; // default value
      loadParameter<int>(ar, "hasAngles", hasAngles_, false); // optional
      if (hasAngles_) {
         if (nAtom_ < 3) {
            UTIL_THROW("Error: Cannot have angles with nAtom < 3");
         }
         nAngle_ = nAtom_ - 2;
         loadParameter<int>(ar, "angleType", angleType_);
      } else {
         nAngle_ = 0;
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      hasDihedrals_ = 0; // default value
      loadParameter<int>(ar, "hasDihedrals", hasDihedrals_, false); // optional
      if (hasDihedrals_) {
         if (nAtom_ < 4) {
            UTIL_THROW("Error: Cannot have dihedrals with nAtom < 4");
         }
         nDihedral_ = nAtom_ - 3;
         loadParameter<int>(ar, "dihedralType", dihedralType_);
      } else {
         nDihedral_ = 0;
      }
      #endif

      buildLinear();
   }

   /*
   * Save internal state to an archive.
   */
   void Multiblock::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar & blockLengths_;
      ar & atomTypes_;
      ar & blockBegin_;
      ar << bondType_;
      #ifdef SIMP_ANGLE
      Parameter::saveOptional(ar, hasAngles_, true);
      if (hasAngles_ && nAngle_ > 0) {
         ar << angleType_;
      } 
      #endif
      #ifdef SIMP_DIHEDRAL
      Parameter::saveOptional(ar, hasDihedrals_, true);
      if (hasDihedrals_ && nDihedral_ > 0) {
         ar << dihedralType_;
      } 
      #endif
   }

   /* 
   * Return atomTypes_ for every atom.
   */
   int Multiblock::calculateAtomTypeId(int index) const
   {
      int type = -1;
      for (int i = 0; i < nBlock_; ++i) {
         if (index >= blockBegin_[i]) {
            type = atomTypes_[i];
         } else {
            break;
         }
      }
      assert(type >= 0);
      return type;
   }

   /* 
   * Return bondType_ for every bond.
   */
   int Multiblock::calculateBondTypeId(int molBondId) const
   { return bondType_; }

   #ifdef SIMP_ANGLE
   /* 
   * Return angleType_ for every angle.
   */
   int Multiblock::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /* 
   * Return dihedralType_ for every dihedral.
   */
   int Multiblock::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif

} 
