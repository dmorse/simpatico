#ifndef MCMD_HOMORING_H
#define MCMD_HOMORING_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Ring.h"

namespace McMd
{

   using namespace Util;

   /**
   * A HomoRing species of loop molecules.
   *
   * \ingroup McMd_Species_Module
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

      #ifdef INTER_ANGLE
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

      #ifdef INTER_DIHEDRAL
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
