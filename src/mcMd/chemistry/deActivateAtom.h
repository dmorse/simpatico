#ifndef MCMD_DE_ACTIVATE_ATOM_H
#define MCMD_DE_ACTIVATE_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Atom;

   /**
   * Temporarily de-activate an atom and all associated groups.
   *
   * \param atom Atom object to be de-activated
   */
   void deActivateAtom(Atom& atom);

   /**
   * Re-activate a temporarily de-activated atom.
   *
   * \param atom Atom object to be re-activated
   */
   void reActivateAtom(Atom& atom);

} 
#endif
