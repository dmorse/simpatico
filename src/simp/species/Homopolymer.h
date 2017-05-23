#ifndef SIMP_HOMOPOLYMER_H
#define SIMP_HOMOPOLYMER_H

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
   * A Homopolymer species of chain molecules.
   *
   * A hompolymer is a chemically homogeneous linear chain of nAtom
   * atoms. All atoms are of the same atom type, with all bonds of 
   * the same bond type. A Homopolymer may optionally have angle
   * potentials for each sequences of three susbseuqent atoms, in
   * which case all angle potentials must be of the same angle type.
   * A Homopolymer may also optionally have dihedral potentials among
   * all sequences of four subsequent atoms, which (if present) must 
   * all be of the same dihedral type.
   *
   * \sa \ref simp_species_Homopolymer_page "parameter file format"
   * 
   * \ingroup Simp_Species_Module
   */
   class Homopolymer : public Linear
   {
   
   public:
  
      /* 
      * Default constructor.
      */
      Homopolymer();
  
      /* 
      * Destructor.
      */
      virtual ~Homopolymer()
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
      int atomType_;
   
      /**
      * Bond type id for every bond of this species.
      */
      int bondType_;
  
      #ifdef SIMP_ANGLE 
      /**
      * Angle type id for every angle of this species.
      */
      int angleType_;
      #endif
   
      #ifdef SIMP_DIHEDRAL
      /**
      * Dihedral type id for every dihedral of this species.
      */
      int dihedralType_;
      #endif
   
      /**
      * Read nAtom_ and the chain type.
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
      * Return the same type for any particle in any chain.
      *
      * \param index atom index, in range 0,...,nAtom_ - 1
      * \return atom type index
      */
      virtual int calculateAtomTypeId(int index) const;
   
      /**
      * Return same bond type for any bond in any chain.
      *
      * \param index bond index, in range 0,...,nBond_ - 1
      * \return bond type index
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef SIMP_ANGLE
      /**
      * Return same angle type for any angle in any chain.
      *
      * \param index  angle index, in range 0,...,nAngle_ - 1
      * \return  angle type index
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Return same dihedral type for any dihedral in any chain.
      *
      * \param index  angle index, in range 0,...,nDihedral_ - 1
      * \return  dihedral type index
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif

   };
   
} 
#endif
