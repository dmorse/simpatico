#ifndef MCMD_HOMOPOLYMER_CPP
#define MCMD_HOMOPOLYMER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
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
      #ifdef INTER_ANGLE
      , angleType_(NullIndex)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralType_(NullIndex)
      #endif
   { setClassName("Homopolymer"); } 
   
   /* 
   * Read nAtom and type.
   */
   void Homopolymer::readSpeciesParam(std::istream &in)
   {
      read<int>(in,"nAtom", nAtom_);
      read<int>(in,"atomType", atomType_);
      nBond_  = nAtom_ - 1;
      read<int>(in,"bondType", bondType_);

      #if INTER_ANGLE
      hasAngles_ = 0;
      read<int>(in, "hasAngles", hasAngles_, false); // optional
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

      #if INTER_DIHEDRAL
      hasDihedrals_ = 0;
      read<int>(in, "hasDihedrals", hasDihedrals_, false); // optional
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
   void Homopolymer::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar,"nAtom", nAtom_);
      loadParameter<int>(ar,"atomType", atomType_);
      nBond_  = nAtom_ - 1;
      loadParameter<int>(ar,"bondType", bondType_);
      #ifdef INTER_ANGLE
      hasAngles_ = 0;
      loadParameter<int>(ar,"hasAngles", hasAngles_, false);
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
      hasDihedrals_ = 0;
      loadParameter<int>(ar,"hasDihedrals", hasDihedrals_, false);
      if (hasDihedrals_) {
         if (nAtom_ > 3) {
            nDihedral_ = nAtom_ - 3;
         } else {
            nDihedral_ = 0;
         }
         if (nDihedral_ > 0) {
            loadParameter<int>(ar, "dihedralType", dihedralType_);
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
   void Homopolymer::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar << nAtom_;
      ar << atomType_;
      ar << bondType_;
      #ifdef INTER_ANGLE
      Parameter::saveOptional(ar, hasAngles_, hasAngles_);
      if (hasAngles_ && nAngle_ > 0) {
         ar << angleType_;
      } 
      #endif
      #ifdef INTER_DIHEDRAL
      Parameter::saveOptional(ar, hasDihedrals_, hasDihedrals_);
      if (hasDihedrals_ && nDihedral_ > 0) {
         ar << dihedralType_;
      } 
      #endif
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

   #ifdef INTER_ANGLE
   /* 
   * Return angleType_ for every angle.
   */
   int Homopolymer::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /* 
   * Return dihedralType_ for every dihedral.
   */
   int Homopolymer::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif
 
} 
#endif
