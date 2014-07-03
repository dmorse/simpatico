#ifndef DDMD_ATOM_CONTEXT_H
#define DDMD_ATOM_CONTEXT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#ifdef DDMD_MOLECULES

namespace DdMd{
   /**
   * Descriptor for context of an Atom within a molecule and species.
   *
   * The information in the AtomContext struct is not used by any 
   * of of the essential algorithms of ddSim, but may be used in the
   * implementation of some analyzers or modifiers.
   */
   struct AtomContext{

      /**
      * Index of species of molecule containing atom.
      */
      int speciesId;

      /**
      * Index of molecule within its molecular species.
      */
      int moleculeId;

      /**
      * Index of atom within its parent molecule.
      *
      * Note: This takes on values from zero to the number of
      * of atoms within the molecule. Corresponding atoms in
      * identical molecules have the same atomId. This index is
      * thus distinct from the global atom id, which is unique
      * within the system.
      */
      int atomId;

   };
}

#endif
#endif
