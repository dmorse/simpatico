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
   struct AtomContext{
      int speciesId;
      int moleculeId;
      int atomId;
   };
}

#endif

#endif
