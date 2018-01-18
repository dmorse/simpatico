/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HomopolymerSG.h"
#ifdef UTIL_MPI
#include <mcMd/simulation/McMd_mpi.h>
#endif

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor.
   */
   HomopolymerSG::HomopolymerSG()
    : Linear(),
      SpeciesMutator(),
      bondType_(NullIndex)
      #ifdef SIMP_ANGLE
      , angleType_(NullIndex)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralType_(NullIndex)
      #endif
   {
      setMutatorPtr(this);
      setClassName("HomopolymerSG");
   } 
   
   /* 
   * Destructor.
   */
   HomopolymerSG::~HomopolymerSG()
   {}

   /*
   * Read a pair of monomer types and exchange chemical potential.
   */
   void HomopolymerSG::readParameters(std::istream& in)
   {
      Species::readParameters(in);
   }

   /* 
   * Read a pair of monomer types and exchange chemical potential.
   */
   void HomopolymerSG::readSpeciesParam(std::istream& in)
   {
      read<int>(in, "nAtom", nAtom_);
      read<int>(in, "bondType", bondType_);
      nBond_ = nAtom_ - 1;

      #ifdef SIMP_ANGLE
      hasAngles_ = 0;  // Default value
      nAngle_ = 0;
      readOptional<int>(in, "hasAngles", hasAngles_);
      if (hasAngles_) {
         if (nAtom_ < 3) {
            UTIL_THROW("Error: Cannot have angles with nAtom < 3");
         }
         nAngle_ = nAtom_ - 2;
         read<int>(in, "angleType", angleType_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      hasAngles_ = 0;  // Default value
      nDihedral_ = 0;  // Default value
      readOptional<int>(in, "hasDihedrals", hasDihedrals_);
      if (hasDihedrals_) {
         if (nAtom_ < 4) {
            UTIL_THROW("Error: Cannot have angles with nAtom < 4");
         }
         nDihedral_ = nAtom_ - 3;
         read<int>(in, "angleType", angleType_);
      }
      #endif

      buildLinear();

      read<Pair <int> >(in, "typeIds", typeIds_);
      read<double>(in, "weightRatio", weightRatio_);

      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }
   
   /* 
   * Load from Serializable::IArchive.
   */
   void HomopolymerSG::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar,"nAtom", nAtom_);
      nBond_  = nAtom_ - 1;
      loadParameter<int>(ar,"bondType", bondType_);

      #ifdef SIMP_ANGLE
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

      #ifdef SIMP_DIHEDRAL
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

      allocateSpeciesMutator(capacity(), 2);

      loadParameter<Pair <int> >(ar, "typeIds", typeIds_);
      loadParameter<double>(ar, "weightRatio", weightRatio_);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   /*
   * Save internal state to an archive.
   */
   void HomopolymerSG::save(Serializable::OArchive &ar)
   {
      ar << moleculeCapacity_;
      ar << nAtom_;
      ar << bondType_;
      #ifdef SIMP_ANGLE
      Parameter::saveOptional(ar, hasAngles_, hasAngles_);
      if (hasAngles_ && nAngle_ > 0) {
         ar << angleType_;
      } 
      #endif
      #ifdef SIMP_DIHEDRAL
      Parameter::saveOptional(ar, hasDihedrals_, hasDihedrals_);
      if (hasDihedrals_ && nDihedral_ > 0) {
         ar << dihedralType_;
      } 
      #endif
      ar << typeIds_;
      ar << weightRatio_;
   }

   /* 
   * Return NullIndex for every atom.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateAtomTypeId(int index) const
   { return NullIndex; }

   /* 
   * Return 0 for every bond.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateBondTypeId(int index) const
   { return bondType_; }

   #ifdef SIMP_ANGLE
   /* 
   * Return 0 for every angle.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /* 
   * Return 0 for every dihedral.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif
 
   /*
   * Change the type of a specific molecule.
   */
   void HomopolymerSG::setMoleculeState(Molecule& molecule, int stateId)
   {
      int nAtom  = molecule.nAtom();
      for (int i = 0; i < nAtom; ++i) {
         molecule.atom(i).setTypeId(typeIds_[stateId]);
      }
      setMoleculeStateId(molecule, stateId);
   }

} 
