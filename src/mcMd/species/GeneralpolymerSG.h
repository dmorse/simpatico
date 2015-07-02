#ifndef MCMD_GENERALPOLYMER_SG_H
#define MCMD_GENERALPOLYMER_SG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Linear.h"                 // base class
#include "SpeciesMutator.h"         // base class
#include <util/containers/DArray.h> // member template
#include <util/containers/Pair.h>   // member

namespace McMd
{

   using namespace Util;

   /**
   * A Homopolymer with a mutable type, for semigrand ensemble.
   *
   * A HomopolymerSG molecule is a Linear chain in which all atoms
   * are of the same type, as are all bonds, but in which the type
   * index for all atoms can be toggled between two values.
   *
   * \ingroup McMd_Species_Module
   */
   class GeneralpolymerSG : public Linear, public SpeciesMutator
   {
   
   public:
  
      /**
      * Constructor.
      */
      GeneralpolymerSG();
   
      /** 
      * Destructor.
      */
      virtual ~GeneralpolymerSG();
  
      /**
      * Read parameters and initialize.
      */ 
      virtual void readParameters(std::istream& in);

      /**
      * Set the type of all atoms in the molecule.
      *
      * \param molecule reference to molecule of interest.
      * \param stateId  atom type id.
      */
      virtual void setMoleculeState(Molecule& molecule, int stateId);

      /**
      * Generate random molecules
      *
      * \param nMolecule number of molecules to genearte
      * \param exclusionRadius array of exclusion radii for every atom type
      * \param system the System
      * \param bondPotentialPtr the bond potential
      * \param boundary the boundary to generate atoms in
      */
      virtual void generateMolecules(int nMolecule,
         DArray<double> exclusionRadius, System &system,
         BondPotential *bondPotentialPtr, const Boundary &boundary);

   protected:
   
      /**
      * Read nAtom, a pair of atom type ids and weightRatio.
      *
      * \param in input stream
      */
      virtual void readSpeciesParam(std::istream &in);
   
      /**
      * Return the same type for any particle in any chain.
      *
      * \param index atom index, in range 0,...,nAtom_ - 1
      * \return atom type index
      */
      virtual int calculateAtomTypeId(int index) const;
   
      /**
      * Return same bond type for any bond in any chain.
      *
      * \param index bond index, in range 0,...,nAtom_ - 1
      * \return bond type index
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef INTER_ANGLE
      /**
      * Return same angle type for any angle in any chain.
      *
      * \param index  angle index, in range 0, ..., nAngle_ - 1
      * \return       angle type index
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Return same dihedral type for any dihedral in any chain.
      *
      * \param index  dihedral index, in range 0, ..., nDihedral_ - 1
      * \return       dihedral type index
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif


   private:

      // A Pair of atom type indexes.
      // Atoms in a molecule with stateId = i have type typeIds_[i]
      Pair<int> typeIds_;

      /// Ratio of statistical weights for typeId[0]/typeId[1]
      double weightRatio_;
      
      /**
      * Try to place an atom. If successful, recursively call tryPlace
      * again to place next atom
      *
      * \param atomIter iterator pointing to the atom to be placed
      * \param exclusionRadius the exclusion radius
      * \param system the System
      * \param cellList the cell list
      * \param bondPotential the bond potential
      *
      * \returns true if particle could be placed
      */
      bool tryPlaceAtom(Molecule &molecule, int atomId,
         DArray<double> exclusionRadius, System& system, CellList &cellList,
         BondPotential *bondPotentialPtr, const Boundary &boundary);

      int maxPlacementAttempts_ = 500;
      DArray<int> beadIdentities_;
      DArray<int> beadTypeIds1_;
      DArray<int> beadTypeIds2_; 
   };
   
} 
#endif
