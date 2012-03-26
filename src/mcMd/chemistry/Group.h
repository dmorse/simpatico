#ifndef MCMD_GROUP_H
#define MCMD_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"

namespace McMd
{

   using namespace Util;

   /**
   * A group of covalently interacting atoms.
   *
   * A group of atoms that interact via a permanent (covalent) 
   * NAtom Atom potential.  Specializations of Group with NAtom=2,
   * 3, and 4 are used to represent groups that interact via
   * covalent bond, angle and dihedral interaction potentials,
   * respectively. The Bond, Angle, and Dihedral typedefs are
   * aliases for these template specializations.
   *
   * A Group<NAtom> contains an array of pointers to Atoms in the 
   * group, and an integer type Id for the group. Parameters of the
   * corresponding Natom body potential are different for groups
   * with different type Ids. 
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
      Group()
       : typeId_(-1),
         isActive_(true)
      {
         for (int i=0; i < NAtom; ++i) {
            atoms_[i] = 0;
         }
      }
            
      /**
      * Get a specific Atom in the Group.
      *
      * \param i index within group (0,...,N-1)
      */
      Atom& atom(int i)
      { return *atoms_[i]; }
    
      /**
      * Get a const reference to a specific Atom in the Group.
      *
      * \param i index within group (0,...,N-1)
      */
      const Atom& atom(int i) const
      { return *atoms_[i]; }
    
      /**
      * Get the typeId for this covalent group.
      */
      int typeId() const
      { return typeId_; }
   
      /**
      * Is this group active?
      *
      * Atoms and the groups containing them can be de-activated and
      * re-activated by some Monte Carlo moves, such as configuration 
      * bias moves that remove and regrow parts of a molecule.
      */
      bool isActive() const
      { return isActive_; }
   
      /**
      * Add an atom to the group.
      *
      * \param i    index of atom within group.
      * \param atom atom to be added.
      */
      void setAtom(int i, Atom &atom)
      { atoms_[i] = &atom; }
     
      /**
      * Set the group type id for this group.
      *
      * \param typeId value of covalent group typeId
      */
      void setTypeId(int typeId)
      { typeId_ = typeId; }
   
      /**
      * Activate or deactivate the group.
      *
      * \param isActive true (active) or false (inactive)
      */
      void setIsActive(bool isActive)
      { isActive_ = isActive; }
   
   private:
      
      /// Array of pointers to Atoms in this group.
      Atom*  atoms_[NAtom];
   
      /// Integer index for the type of group.
      int    typeId_;
   
      /// If false, the group has been temporarily deactivated.
      bool   isActive_;
   
   };

}
#endif
