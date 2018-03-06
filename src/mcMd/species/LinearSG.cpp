/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearSG.h"
#ifdef UTIL_MPI
#include <mcMd/simulation/McMd_mpi.h>
#endif

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LinearSG::LinearSG()
    : Linear(),
      SpeciesMutator(),
      beadTypeIds0_(),
      beadTypeIds1_()
      #ifdef SIMP_ANGLE
      , angleType_(NullIndex)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralType_(NullIndex)
      #endif
   {
      setMutatorPtr(this);
      setClassName("LinearSG");
   }

   /*
   * Destructor.
   */
   LinearSG::~LinearSG()
   {}

   /*
   * Call general Species::readParameters() .
   */
   void LinearSG::readParameters(std::istream& in)
   {
      Species::readParameters(in);
   }
  

   /*
   * Read atom structure and two sets of atom type ids.
   */
   void LinearSG::readSpeciesParam(std::istream& in)
   {
      read<int>(in,"nAtom", nAtom_);
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
      beadTypeIds0_.allocate(nAtom_);
      beadTypeIds1_.allocate(nAtom_);
      readDArray<int>(in, "identities0", beadTypeIds0_, nAtom_);
      readDArray<int>(in, "identities1", beadTypeIds1_, nAtom_);

      read<double>(in, "weightRatio", weightRatio_);

      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   void LinearSG::loadSpeciesParam(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar,"nAtom", nAtom_);
      loadParameter<int>(ar,"bondType",bondType_);
      nBond_  = nAtom_ - 1;

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

      loadParameter<Pair <int> >(ar, "typeIds", typeIds_);
      
      ar & beadTypeIds0_;
      ar & beadTypeIds1_;

      loadParameter<double>(ar, "weightRatio", weightRatio_);
      allocateSpeciesMutator(capacity(), 2);

      // Set statistical weights
      double sum = weightRatio_ + 1.0;
      mutator().setWeight(0, weightRatio_/sum);
      mutator().setWeight(1, 1.0/sum);
   }

   /*
   * Save internal state to an archive.
   */
   void LinearSG::save(Serializable::OArchive &ar)
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
      ar << beadTypeIds0_;
      ar << beadTypeIds1_;
      ar << weightRatio_;
   }

   /*
   * Return NullIndex for every atom.
   * Set initial typeIds
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAtomTypeId(int index) const
   { return NullIndex; }

   /*
   * Return 0 for every bond.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateBondTypeId(int index) const
   { return bondType_; }

   #ifdef SIMP_ANGLE
   /*
   * Return 0 for every angle.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateAngleTypeId(int index) const
   { return angleType_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Return 0 for every dihedral.
   *
   * Used by Linear::buildLinear().
   */
   int LinearSG::calculateDihedralTypeId(int index) const
   { return dihedralType_; }
   #endif

   /*
   * Change the type of a specific molecule.
   */
   void LinearSG::setMoleculeState(Molecule& molecule, int stateId)
   {
      int nAtom  = molecule.nAtom();
      DArray<int> beadIdentities;
      for (int i = 0; i < nAtom; ++i) {
         if (stateId == 0) {
            beadIdentities = beadTypeIds0_;
         } else {
            beadIdentities = beadTypeIds1_;
         }
         molecule.atom(i).setTypeId(beadIdentities[i]);
      }
      setMoleculeStateId(molecule, stateId);
   }

}
