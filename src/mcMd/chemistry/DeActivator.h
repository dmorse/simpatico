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
   * Class of static methods to de-active and re-active atoms.
   *
   * Static methods of this class can de-activate and re-activate
   * individual atoms or all atoms in a molecule, while guaranteeing
   * the consistency of the inactive atom counter in all associated
   * Group objects. 
   *
   * This class is a friend of Atom and Group<N> classes to allow 
   * access to private functions that modify Atom::Active_ flag 
   * and the Group<N>::nInActive_ inactive atom counter.
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
      * Sets isActive flag to true for all atoms and inactive atom
      * counter to zero for all groups.
      *
      * \param molecule Molecule object to be fully activated
      */
      static void activate(Molecule& molecule);
   
   };

} 
#endif
