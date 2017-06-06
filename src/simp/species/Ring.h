#ifndef SIMP_RING_H
#define SIMP_RING_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*
* Author: Jian Qin
*/

#include "Species.h"

namespace Simp
{

   using namespace Util;

   /**
   * A Species of ring polymers (abstract).
   *
   * A Ring molecule is a loop of nMolParticle atoms with consecutive
   * atom indices, with nMolParticle bonds between neighboring atoms.
   *
   * A Ring molecule may be chemically heterogeneous - different atoms and 
   * bonds in a loop may be of different types. The chemical structure of a 
   * subclass of Ring must be specified by implementing the pure virtual 
   * functions atomTypeId() and bondTypeId().
   *
   * \ingroup Simp_Species_Module
   */
   class Ring : public Species 
   {
   
   public:
   
      // Methods
   
      /**
      * Constructor.
      */
      Ring(); 

      /**
      * Destructor.
      */
      virtual ~Ring();

   protected:

      /**
      * Return the atom type id for a specific atom.
      *
      * \param index atom index, in the range 0, ... , nAtom_-1
      * \return type of the specified bond
      */
      virtual int calculateAtomTypeId(int index) const = 0;

      /**
      * Return the bond type id for a specific bond.
      *
      * Bond i of a ring molecule connects atom i and atom i+1.
      *
      * \param index bond index, in the range 0, ... , nBond_-1
      * \return type of the specified bond
      */
      virtual int calculateBondTypeId(int index) const = 0;

      #ifdef SIMP_ANGLE
      /**
      * Return the angle type id for a specific angle.
      *
      * Angle i of a ring molecule connects atoms: (i)-(i+1)-(i+2), where
      * atom indices mod to nAtom_.
      *
      * \param index local angle index, in the range 0, ... , nAngle_ - 1
      * \return type of the specified bond
      */
      virtual int calculateAngleTypeId(int index) const = 0;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Return the dihedral type id for a specific dihedral.
      *
      * Dihedral i of a ring molecule connects atoms: (i)-(i+1)-(i+2)-(i+3),
      * where atom indices mod to nAtom_.
      *
      * \param index local dihedral index, in the range 0, ... , nDihedral_ - 1
      * \return type of the specified bond
      */
      virtual int calculateDihedralTypeId(int index) const = 0;
      #endif

      /**
      * Build the chemical structure for a ring molecule.
      *
      * This method must be called by readSpeciesParam() after values are 
      * set for nAtom and nBond. It allocates memory, assigns atom types to
      * all the atoms in a generic molecule, and creates SpeciesGroup<2>
      * objects for all of the bonds.
      */
      void buildRing();
   
   };

} 
#endif
