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

      /**
      * Activate or deactivate the group.
      *
      * \param isActive true (active) or false (inactive)
      */
      void setIsActive(bool isActive);

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
      * Atoms and the groups containing them can be de-activated and
      * re-activated by some Monte Carlo moves, such as configuration 
      * bias moves that remove and regrow parts of a molecule.
      */
      bool isActive() const;

   private:

      /// Array of pointers to Atoms in this group.
      Atom* atoms_[NAtom];

      /// Integer index for the type of group.
      int typeId_;

      /// If false, the group has been temporarily deactivated.
      bool isActive_;

   };

   /*
   * Constructor 
   */
   template <int NAtom>
   Group<NAtom>::Group()
    : typeId_(-1),
      isActive_(true)
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

   /*
   * Activate or deactivate the group.
   */
   template <int NAtom>
   inline void Group<NAtom>::setIsActive(bool isActive)
   {  isActive_ = isActive; }

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
   {  return isActive_; }

}
#endif
