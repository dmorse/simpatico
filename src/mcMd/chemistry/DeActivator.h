#ifndef MCMD_DE_ACTIVATOR_H
#define MCMD_DE_ACTIVATOR_H

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
   * This class has only static member functions, and no member 
   * variables. It provides functions that can de-activate and
   * re-activate individual atoms and activate all atoms in a
   * molecule, while guaranteeing the consistency of the count
   * of inactive atoms in associated covalent group objects,
   * i.e., all bonds, angles, and dihedrals that contain that
   * atom.
   *
   * This set of functions is defined as class to allow a
   * friend class declaration: Class DeActivator is declared 
   * a friend of Atom and Group<N> classes to allow access to 
   * private member functions that modify the Atom::isActive_ 
   * bool flag and the Group<N>::nInActive_ int counter. By 
   * allowing these functions to be called only via members
   * of DeActivator, which call all required functions when 
   * the status of an atom changes, we guarantee that Atom 
   * and Group<N> objects are in consistent states.
   */
   class DeActivator
   {
 
   public: 

      /**
      * Temporarily de-activate one atom and all associated groups.
      *
      * \param atom Atom object to be de-activated
      */
      static void deActivate(Atom& atom);
   
      /**
      * Re-activate a temporarily de-activated atom and update groups.
      *
      * \param atom Atom object to be re-activated
      */
      static void reActivate(Atom& atom);

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
