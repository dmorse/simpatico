/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HomopolymerSG.h"
#ifdef UTIL_MPI
#include <mcMd/simulation/mcMd_mpi.h>
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
      SpeciesMutator()
   {
      setMutatorPtr(this);
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
      read<int>(in,"nAtom", nAtom_);
      nBond_ = nAtom_ - 1;
      #ifdef SIMP_ANGLE
      nAngle_ = nAtom_ - 2;
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nAtom_ > 3)
         nDihedral_ = nAtom_ - 3;
      else
         nDihedral_ = 0;
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
   { return 0; }

   #ifdef SIMP_ANGLE
   /* 
   * Return 0 for every angle.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateAngleTypeId(int index) const
   { return 0; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /* 
   * Return 0 for every dihedral.
   *
   * Used by Linear::buildLinear().
   */
   int HomopolymerSG::calculateDihedralTypeId(int index) const
   { return 0; }
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
