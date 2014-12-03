#ifndef MCMD_ACTIVATE_H
#define MCMD_ACTIVATE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Atom;
   class Molecule;

   /**
   * Static member functions to de-active and re-active atoms.
   *
   * This class has only 3 static member functions, and no member 
   * variables. It provides functions that can de-activate and
   * re-activate individual atoms or activate all atoms in a
   * molecule. The implementation  consistently modifies the count 
   * of inactive atoms in all associated covalent group objects,
   * i.e., in all bonds, angle, and dihedral groups that contain
   * an atom.
   *
   * This set of functions is defined as class (rather than 
   * simply a set of functions) to allow use of friend class 
   * declarations: Class Activate is declared a friend of the
   * Atom class and the Group<N> class template. This allows 
   * access to private member functions that can modify the 
   * Atom::isActive_ bool flag and the Group<N>::nInActive_ 
   * integer counter. By allowing these private functions to 
   * be called only by members of Activate, we can guarantee 
   * that associated Atom and Group<N> objects maintain 
   * consistent states.
   * 
   * \ingroup McMd_Chemistry_Module
   */
   class Activate
   {
 
   public: 

      /**
      * Temporarily de-activate one atom and update associated groups.
      *
      * \param atom Atom object to be de-activated
      */
      static void deactivate(Atom& atom);
   
      /**
      * Re-activate a temporarily de-activated atom and update groups.
      *
      * \param atom Atom object to be re-activated
      */
      static void reactivate(Atom& atom);

      /**
      * Activate all atoms and groups in one molecule.
      *
      * Marks all atoms as active and sets the number of inactive 
      * atoms to zero for all groups in the molecule.
      *
      * \param molecule Molecule object to be fully activated
      */
      static void activate(Molecule& molecule);
   
   };

} 
#endif
