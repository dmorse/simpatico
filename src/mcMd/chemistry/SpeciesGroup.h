#ifndef MCMD_SPECIES_GROUP_H
#define MCMD_SPECIES_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <iostream>

namespace McMd
{

   // Forward declarations of templates required for friend declarations.

   template <int NAtom> class SpeciesGroup;

   template <int NAtom> std::istream& 
   operator >> (std::istream& in, SpeciesGroup<NAtom> &speciesGroup);

   template <int NAtom> std::ostream& 
   operator << (std::ostream& out, const SpeciesGroup<NAtom>& speciesGroup);

   template <int NAtom> std::ostream& 
   operator << (std::ostream& out, const SpeciesGroup<NAtom>& speciesGroup);

   /**
   * Serialize one SpeciesGroup<NAtom>.
   *
   * Default implementation calls serialize member function of data object.
   * Can be overridden by any explicit specialization.
   *
   * \param ar            archive object
   * \param speciesGroup  object to be serialized
   * \param version       archive version id
   */
   template <class Archive, int NAtom>
   void serialize(Archive& ar, SpeciesGroup<NAtom>& speciesGroup, 
                  const unsigned int version);

   /**
   * A Group of covalently interacting atoms within any molecule of one Species.
   *
   * Class SpeciesGroup is used in the implementation of class Species to define 
   * the covalent chemical structure for any molecule of a particular species. 
   * The chemical structure of a Species is defined by a set of SpeciesGroup<2> 
   * (bond), SpeciesGroup<3> (angle), and SpeciesGroup<4> (dihedral) objects.
   * Each SpeciesGroup<NAtom> object contains an array of NAtom local integer
   * atom indices, which identities a group of atoms that are involved in a
   * an NAtom-body covalent interaction. Each local atom index i must lie in 
   * the range 0 <= i < nAtom. A local atom index identifies the position of
   * of an Atom within a generic molecule of the associated species. Each 
   * SpeciesGroup also has an integer group type index to distinguish different 
   * types of bonds, angles, and dihedrals, with different potential parameters.
   *
   * The difference between a SpeciesGroup and a Group is that a Group contains
   * an array of pointers to specific Atom objects that all belong to a specific 
   * Molecule. Each SpeciesGroup is thus a prototype for constructing one Group 
   * for each Molecule within a Species.
   *
   * \ingroup McMd_Chemistry_Module
   */
   template <int NAtom>
   class SpeciesGroup
   {

   public: 
     
      /**
      * Null (unknown) value for non-negative atom and group type indices.
      */
      static const int NullIndex = -1;

      /**
      * Constructor.
      */
      SpeciesGroup() 
       : typeId_(NullIndex)
      {
         for (int i=0; i < NAtom; ++i) {
            atomIds_[i] = NullIndex;
         }
      }
   
      /**
      * Set index for one atom in the group.
      *
      * \param i   index of atom within group (0, .., NAtom - 1)
      * \param Id  index of atom to be added.
      */
      void setAtomId(int i, int Id)
      { atomIds_[i] = Id; }
     
      /**
      * Set the group type id for this group.
      *
      * \param typeId value of covalent group typeId
      */
      void setTypeId(int typeId)
      { typeId_ = typeId; }

      /**
      * Get the local id for a specific Atom.
      *
      * \param i index within group (0, ..., NAtom - 1)
      */
      int atomId(int i) const
      { return atomIds_[i]; }
    
      /**
      * Get the type index for this covalent group.
      */
      int typeId() const
      { return typeId_; }
   
      #ifdef UTIL_MPI

      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();

      #endif

   private:
      
      /// Array of atom Ids in this group.
      int    atomIds_[NAtom];
   
      /// Integer index for the type of group.
      int    typeId_;
   
   //friends:

      friend
      std::istream& 
      operator >> <> (std::istream& in, SpeciesGroup<NAtom>& speciesGroup);

      friend 
      std::ostream& 
      operator << <> (std::ostream& out, const SpeciesGroup<NAtom>& speciesGroup);

      template <class Archive> friend
      void serialize(Archive& ar, SpeciesGroup<NAtom>& speciesGroup, 
                     const unsigned int version);

      template <class Archive> friend
      void serialize(Archive& ar, SpeciesGroup<NAtom>& speciesGroup, 
                     const unsigned int version);

   };

   // Inserter and extractor operator declarations

   /**
   * istream extractor for a SpeciesGroup.
   *
   * The format is as follows:
   *
   *      atomIds[0] atomIds[1] ...  atomId[NAtom-1]  typeId
   *
   * \param in            input stream
   * \param speciesGroup  SpeciesGroup to be read from stream
   * \return modified input stream
   */
   template <int NAtom>
   std::istream& 
   operator >> (std::istream& in, SpeciesGroup<NAtom> &speciesGroup);

   /**
   * ostream inserter for a SpeciesGroup.
   *
   * The format is as follows, output one one line with no line break:
   *
   *      atomIds[0] atomIds[1] ...  atomId[NAtom-1] typeId
   *
   * \param  out           output stream
   * \param  speciesGroup  SpeciesGroup to be written to stream
   * \return modified output stream
   */
   template <int NAtom>
   std::ostream& 
   operator << (std::ostream& out, const SpeciesGroup<NAtom> &speciesGroup);

   template <class Archive>
   void serialize(Archive& ar, SpeciesGroup<2>& speciesGroup, const unsigned int version)
   {
      for (int i = 0; i < 2; ++i) {
         ar & speciesGroup.atomIds_[i];
      }
      ar & speciesGroup.typeId_;
   }

   template <class Archive>
   void serialize(Archive& ar, SpeciesGroup<3>& speciesGroup, const unsigned int version)
   {
      for (int i = 0; i < 3; ++i) {
         ar & speciesGroup.atomIds_[i];
      }
      ar & speciesGroup.typeId_;
   }

   template <class Archive>
   void serialize(Archive& ar, SpeciesGroup<4>& speciesGroup, const unsigned int version)
   {
      for (int i = 0; i < 4; ++i) {
         ar & speciesGroup.atomIds_[i];
      }
      ar & speciesGroup.typeId_;
   }

}
#endif
