#ifndef MCMD_GROUP_H
#define MCMD_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * A sequence of NAtom covalently interacting atoms.
   *
   * An ordered group of atoms that interact via a permanent 
   * (covalent) NAtom atom potential energy.  A Group<NAtom> 
   * object has an array of pointers to Atoms in the group.
   * Specializations of Group<NAtom> with NAtom=2, 3, and 4 
   * are used to represent covalent bond, angle and dihedral 
   * interaction groups, respectively. The Bond, Angle, and 
   * Dihedral typedefs are aliases for these template 
   * specializations.
   *
   * A Group<NAtom> object also has an integer type id. 
   * Different parameters may be used in the corresponding 
   * covalent potential for groups with different type Ids. 
   *
   * \ingroup McMd_Chemistry_Module
   */
   template <int NAtom>
   class Group
   {

   public: 

      /**
      * Constructor 
      */
      Group();

      /// \name Initialization
      //@{

      /**
      * Add an atom to this group.
      *
      * \param i    index of atom within group.
      * \param atom atom to be added.
      */
      void setAtom(int i, Atom &atom);

      /**
      * Set the group type id for this group.
      *
      * \param typeId value of covalent group typeId
      */
      void setTypeId(int typeId);

      //@}
      /// \name Accessors
      //@{

      /**
      * Get a specific Atom in the Group by reference.
      *
      * \param i index within group (0,...,N-1)
      */
      Atom& atom(int i);

      /**
      * Get a specific Atom in the Group by const reference.
      *
      * \param i index within group (0,...,N-1)
      */
      const Atom& atom(int i) const;

      /**
      * Get the typeId for this covalent group.
      */
      int typeId() const;

      /**
      * Is this group active?
      *
      * A group is active iff all associated atoms are active. Atoms
      * can be temporarily de-activated during some Monte Carlo moves
      * that remove and regrow parts of a molecule. 
      */
      bool isActive() const;

      /**
      * Return number of inactive atoms.
      */
      int nInActive() const;

      /**
      * Check consistency of number of inactive atoms.
      *
      * Returns true if consistent, or throws Exception.
      */
      bool checkInactive() const;

      //@}

   private:

      /// Array of pointers to Atoms in this group.
      Atom* atoms_[NAtom];

      /// Integer index for the type of group.
      int typeId_;

      /// Number of inactive atoms in group.
      int nInActive_;

      // Private member functions, accessible by class Activate

      /**
      * Activate this group (set number of inactive atoms to zero).
      */
      void activate();

      /**
      * Increment the number of inactive atoms.
      */
      void incrementInactive();

      /**
      * Decrement the number of inactive atoms.
      */
      void decrementInactive();

   // friends:

      friend class Activate;

   };

   /*
   * Constructor 
   */
   template <int NAtom>
   Group<NAtom>::Group()
    : typeId_(-1),
      nInActive_(0)
   {
      for (int i=0; i < NAtom; ++i) {
         atoms_[i] = 0;
      }
   }

   /*
   * Add an atom to the group.
   */
   template <int NAtom>
   inline void Group<NAtom>::setAtom(int i, Atom &atom)
   {  atoms_[i] = &atom; }

   /*
   * Set the group type id for this group.
   */
   template <int NAtom>
   inline void Group<NAtom>::setTypeId(int typeId)
   {  typeId_ = typeId; }

   // Accessors

   /*
   * Get a specific Atom in the Group.
   */
   template <int NAtom>
   inline Atom& Group<NAtom>::atom(int i)
   {
      assert(atoms_[i]);  
      return *atoms_[i]; 
   }

   /*
   * Get a const reference to a specific Atom in the Group.
   */
   template <int NAtom>
   inline const Atom& Group<NAtom>::atom(int i) const
   {  
      assert(atoms_[i]);  
      return *atoms_[i]; 
   }

   /*
   * Get the typeId for this covalent group.
   */
   template <int NAtom>
   inline int Group<NAtom>::typeId() const
   {  return typeId_; }

   /*
   * Is this group active?
   */
   template <int NAtom>
   inline bool Group<NAtom>::isActive() const
   {  return (nInActive_ == 0); }

   /*
   * Return number of inactive atoms.
   */
   template <int NAtom>
   inline int Group<NAtom>::nInActive() const
   {  return nInActive_; }

   /*
   * Check consistency of nInActive counter.
   */
   template <int NAtom>
   bool Group<NAtom>::checkInactive() const
   {
      int counter = 0;
      for (int i=0; i < NAtom;  ++i) {
         assert(atoms_[i]);  
         if (!atoms_[i]->isActive()) ++counter;
      }
      if (counter != nInActive_) {
         UTIL_THROW("Inconsistent number of inactive atoms");
      }
      return true;
   }

   // Private functions that modify nInActive_

   /*
   * Activate the group (set nInActive_ = 0)   
   */
   template <int NAtom>
   inline void Group<NAtom>::activate()
   {  nInActive_ = 0; }

   /*
   * Increment the number of inactive atoms.
   */
   template <int NAtom>
   inline void Group<NAtom>::incrementInactive()
   {  
      assert(nInActive_ < NAtom);  
      ++nInActive_; 
   }

   /*
   * Decrement the number of inactive atoms.
   */
   template <int NAtom>
   inline void Group<NAtom>::decrementInactive()
   {
      assert(nInActive_ > 0);
      --nInActive_;
   }

}
#endif
