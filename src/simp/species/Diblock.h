#ifndef SIMP_DIBLOCK_H
#define SIMP_DIBLOCK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Linear.h"

namespace Simp
{

   using namespace Util;

   /**
   * A linear diblock copolymer chain.
   *
   * In this implementation, all of the bonds have the same type Id.
   *
   * \sa \ref simp_species_Diblock_page "parmameter file format"
   *
   * \ingroup Simp_Species_Module
   */
   class Diblock : public Linear
   {
   
   public:
   
      /// Default constructor.
      Diblock();
   
      /// Destructor.
      virtual ~Diblock()
      {}
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
   protected:
   
      /// Number of particles in each block.
      int blockLengths_[2];  
   
      /// Particle type ids for each block
      int atomTypes_[2];         
 
      int bondType_;
 
      #ifdef SIMP_ANGLE 
      int angleType_;
      #endif

      #ifdef SIMP_DIHEDRAL
      int dihedralType_;
      #endif

      /**
      * Read block lengths length and and types.
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
      * Return the atom type for a specific atom
      *
      * \param index index of atom within the molecule
      * \return atom type for the specified atom
      */
      virtual int calculateAtomTypeId(int index) const;
   
      /**
      * Return same bond type for any bond in any chain.
      *
      * \param index index of bond within the molecule
      * \return bond type index for the specified bond
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef SIMP_ANGLE
      /**
      * Return same angle type for any angle in any chain.
      *
      * \param index  angle index, in range 0, ..., nAngle_ - 1
      * \return bond  type index
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Return same dihedral type for any dihedral in any chain.
      *
      * \param index  dihedral index, in range 0, ..., nDihedral_ - 1
      * \return  dihedral type index
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif
 
   };
 
} 
#endif
