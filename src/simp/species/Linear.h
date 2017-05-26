#ifndef SIMP_LINEAR_H
#define SIMP_LINEAR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"

namespace Simp
{

   using namespace Util;

   class System;
   class CellList;
   class BondPotential;
   class Molecule;

   /**
   * A Species of linear polymers (abstract).
   *
   * A linear molecule is a chain of nAtom atoms with consecutive atom
   * indices, with nAtom - 1 bonds between neighboring atoms.
   *
   * A Linear molecule may be chemically heterogeneous - different atoms and 
   * bonds in a chain may be of different types. The chemical structure of a 
   * subclass of Linear must be specified by implementing the pure virtual 
   * functions atomTypeId() and bondTypeId().
   *
   * \ingroup Simp_Species_Module
   */
   class Linear : public Species 
   {
   
   public:
   
      /**
      * Constructor.
      */
      Linear(); 

      /**
      * Destructor.
      */
      virtual ~Linear();

   protected:

      #ifdef SIMP_ANGLE 
      /**
      * Does this chain have angle potentials (0 = false, 1 = true).
      */
      int hasAngles_;
      #endif
   
      #ifdef SIMP_DIHEDRAL
      /**
      * Does this chain have dihedral potentials (0 = false, 1 = true).
      */
      int hasDihedrals_;
      #endif
   
      /**
      * Return the atom type id for a specific atom.
      *
      * Implementations of this function should return a type Id for
      * a specified atom that this is calculated from the information 
      * read by readParameters().  It is used in Linear::buildLinear() to 
      * build a chain.
      * 
      * \param index local atom index, in the range 0, ... , nAtom_-1
      * \return type of the specified bond
      */
      virtual int calculateAtomTypeId(int index) const = 0;

      /**
      * Return the bond type id for a specific bond.
      *
      * Bond i of a linear molecule connects atom i and atom i + 1.
      *
      * \param index local bond index, in the range 0, ... , nBond_ - 1
      * \return type of the specified bond
      */
      virtual int calculateBondTypeId(int index) const = 0;

      #ifdef SIMP_ANGLE
      /**
      * Return the angle type id for a specific angle.
      *
      * Angle i of a linear molecule connects atoms: i-(i+1)-(i+2).
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
      * Dihedral i of a ring molecule connects atoms: (i)-(i+1)-(i+2)-(i+3).
      *
      * \param index local dihedral index, in the range 0, ... , nDihedral_ - 1
      * \return type of the specified bond
      */
      virtual int calculateDihedralTypeId(int index) const = 0;
      #endif

      /**
      * Build the chemical structure for a linear molecule.
      *
      * This method must be called by readSpeciesParam() after values are
      * set for nAtom and nBond. It allocates memory, assigns atom types to
      * all the atoms in a generic molecule, and creates SpeciesGroup<2>
      * objects for all of the bonds.
      */
      void buildLinear();

   };

} 
#endif
