#ifndef SIMP_HOMO_RING_H
#define SIMP_HOMO_RING_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*
* Author: Jian Qin
*/

#include "Ring.h"

namespace Simp
{

   using namespace Util;

   /**
   * A species of homogeneous ring molecules.
   *
   * A HomoRing is a chemically homogeneous ring of atoms, in which
   * all atoms are of the same type.
   *
   * \sa \ref simp_species_HomoRing_page "parameter file format"
   * \ingroup Simp_Species_Module
   */
   class HomoRing : public Ring
   {
   
   public:
   
      /// Default constructor.
      HomoRing();
   
      /// Destructor.
      virtual ~HomoRing()
      {}
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
   protected:
   
      /**
      * Particle type id for every  particle of every molecule of this species.
      */
      int type_;
   
      /**
      * Read nAtom_ and the loop type.
      *
      * \param in input stream
      */
      virtual void readSpeciesParam(std::istream &in);
   
      /**
      * Load species structure from an Archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

      /**
      * Return the same type for any particle in any loop.
      *
      * \param index  atom index, in range 0, ..., nAtom_ - 1
      * \return atom  type index
      */
      virtual int calculateAtomTypeId(int index) const;
   
      /**
      * Return same bond type for any bond in any loop.
      *
      * \param index  bond index, in range 0, ..., nBond_ - 1
      * \return bond  type index
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef SIMP_ANGLE
      /**
      * Return same angle type for any homogeneous ring.
      *
      * Angle i of a ring molecule connects atoms: (i)-(i+1)-(i+2), where
      * atom indices mod to nAtom_.
      *
      * \param index local angle index, in the range 0, ... , nAngle_ - 1
      * \return type of the specified angle
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Return same dihedral type for any homogeneous ring.
      *
      * Angle i of a ring molecule connects atoms: (i)-(i+1)-(i+2)-(i+3), where
      * atom indices mod to nAtom_.
      *
      * \param index local dihedral index, in the range 0, ... , nDihedral_ - 1
      * \return type of the specified dihedral
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif

   };
   
} 
#endif
