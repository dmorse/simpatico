#ifndef DDMD_GROUP_H
#define DDMD_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"

namespace DdMd
{

   using namespace Util;

   // Forward declarations required for friend declarations

   template <int N> class Group;

   template <int N>
   std::istream& operator >> (std::istream& in, Group<N> &group);

   template <int N>
   std::ostream& operator << (std::ostream& out, const Group<N> &group);

   template <class Archive, int N>
   void serialize(Archive& ar, Group<N>& group, const unsigned int version);

   /**
   * A group of covalently interacting atoms.
   *
   * A group of atoms that interact via a permanent (covalent) 
   * N Atom potential.  Specializations of Group with N=2, 3,
   * 4 are used to represent groups that interact via covalent
   * bond, angle and torsion interaction potentials, respectively.
   *
   * A Group<N> contains both an array of integer ids for Atoms 
   * in the group, and an array of pointers to these atoms. It
   * also has an integer type Id for the group and a global id. 
   *
   * \ingroup DdMd_Chemistry_Module
   */
   template <int N>
   class Group
   {

   public: 

      /**
      * Constructor 
      */
      Group()
       : typeId_(-1),
         id_(-1),
         nPtr_(0)
      {
         for (int i=0; i < N; ++i) {
            atomIds_[i]  = -1;
            atomPtrs_[i] = 0;
         }
         plan_.clearFlags();
      }
            
      /**
      * Set group to empty initial state.
      */
      void clear()
      {
         typeId_ = -1;
         id_ = -1; 
         nPtr_ = 0;
         for (int i=0; i < N; ++i) {
            atomIds_[i]  = -1;
            atomPtrs_[i] = 0;
         }
         plan_.clearFlags();
      }
   
      /**
      * Set the global id for this group.
      *
      * \param id value of global id (tag) for this group.
      */
      void setId(int id)
      {  id_ = id; }
   
      /**
      * Set the group type id for this group.
      *
      * \param typeId value of covalent group typeId
      */
      void setTypeId(int typeId)
      {  typeId_ = typeId; }
   
      /**
      * Set the id for a specific atom.
      *
      * \param i      index within group (0,..., N-1)
      * \param atomId global atom identifier (atom tag)
      */
      void setAtomId(int i, int atomId) 
      {  atomIds_[i] = atomId; }
    
      /**
      * Set the pointer to a specific atom.
      *
      * \param i       index of atom within group
      * \param atomPtr pointer to Atom to be added
      */
      void setAtomPtr(int i, Atom* atomPtr)
      {
         if (atomPtr == 0) {
            UTIL_THROW("Attempt to set null pointer");
         }  
         if (atomPtrs_[i] == 0) {
            ++nPtr_;
         }
         atomPtrs_[i] = atomPtr; 
      }
    
      /**
      * Clear the pointer to a specific atom.
      *
      * \param i index of atom within group.
      */
      void clearAtomPtr(int i)
      {
         if (atomPtrs_[i] != 0) {
            --nPtr_;
            atomPtrs_[i] = 0; 
         }
      }
  
      /**
      * Get communication plan by reference.
      */
      Plan& plan()
      {  return plan_; }

      // Accessors
 
      /**
      * Get the global id for this group.
      */
      int id() const
      {  return id_; }
   
      /**
      * Get the typeId for this group.
      */
      int typeId() const
      {  return typeId_; }
   
      /**
      * Get the id for a specific atom in the Group.
      *
      * \param i index within group (0,..., N-1)
      */
      int atomId(int i) const
      {  return atomIds_[i]; }
    
      /**
      * Get a pointer to a specific Atom.
      *
      * \param i index within group (0,...,N-1)
      */
      Atom* atomPtr(int i) const
      {  return atomPtrs_[i]; }
    
      /**
      * Return the number of non-null atom pointers in this group.
      */
      int nPtr() const 
      {  return nPtr_; }
   
      /**
      * Get communication plan (const reference).
      */
      const Plan& plan() const
      {  return plan_; }

   private:
      
      /// Array of pointers to Atoms in this group.
      Atom*  atomPtrs_[N];
   
      /// Array of integer ids of atoms in this group.
      int  atomIds_[N];
   
      /// Integer index for the type of group.
      int  typeId_;
   
      /// Global id for this group.
      int  id_;

      /// Number of non-null atom pointers in this Group.
      int  nPtr_;

      // Communication plan.
      Plan plan_;

   //friends:

      friend std::istream& operator >> <> (std::istream& in, Group<N> &group);

      friend std::ostream& operator << <> (std::ostream& out, const Group<N> &group);

      template <class Archive> friend
      void serialize(Archive& ar, Group<N>& group, const unsigned int version);


   };

   // Friend declarations

   /**
   * istream extractor (>>) for a Group.
   *
   * Format:
   *
   * \param in        input stream
   * \param group  Group to be read from stream
   * \return modified input stream
   */
   template <int N>
   std::istream& operator>>(std::istream& in, Group<N> &group);

   /**
   * ostream inserter (<<) for a Group.
   *
   * Format, one one line with no line break:
   *
   * \param  out   output stream
   * \param  group Group to be written to stream
   * \return modified output stream
   */
   template <int N>
   std::ostream& operator<<(std::ostream& out, const Group<N> &group);

   /**
   * Serialize one Group<N>.
   *
   * Default implementation calls serialize method of data object.
   * Can be overridden by any explicit specialization.
   *
   * \param ar            archive object
   * \param group  object to be serialized
   * \param version       archive version id
   */
   template <class Archive, int N>
   void serialize(Archive& ar, Group<N>& group, 
                  const unsigned int version);

   // Friend implementations.

   /* 
   * Input a Goup<N> from an istream, without line breaks.
   */
   template <int N>
   std::istream& operator>>(std::istream& in, Group<N> &group)
   {
      in >> group.id_;
      in >> group.typeId_;
      for (int i=0; i < N; ++i) {
         in >> group.atomIds_[i];
      }
      return in;
   }
   
   /* 
   * Output a Group to an ostream, without line breaks.
   */
   template <int N>
   std::ostream& operator<<(std::ostream& out, const Group<N> &group) 
   {
      out.width(10);
      out << group.id_;
      out.width(10);
      out << group.typeId_;
      for (int i = 0; i < N; ++i) {
         out.width(10);
         out << group.atomIds_[i];
      }
      return out;
   }

   template <class Archive>
   void serialize(Archive& ar, Group<2>& group, const unsigned int version)
   {
      ar & group.id_;
      ar & group.typeId_;
      for (int i = 0; i < 2; ++i) {
         ar & group.atomIds_[i];
      }
   }

   template <class Archive>
   void serialize(Archive& ar, Group<3>& group, const unsigned int version)
   {
      ar & group.id_;
      ar & group.typeId_;
      for (int i = 0; i < 3; ++i) {
         ar & group.atomIds_[i];
      }
   }

   template <class Archive>
   void serialize(Archive& ar, Group<4>& group, const unsigned int version)
   {
      ar & group.id_;
      ar & group.typeId_;
      for (int i = 0; i < 4; ++i) {
         ar & group.atomIds_[i];
      }
   }

} 
#endif
