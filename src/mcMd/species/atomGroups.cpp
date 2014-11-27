/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "atomGroups.h"
#include <mcMd/chemistry/Atom.h>  
#include <mcMd/chemistry/Molecule.h>       

namespace McMd
{

   using namespace Util;

   #ifdef INTER_BOND
   /*
   * Fill an array of pointers to Bonds that contain an Atom.
   */
   void getAtomBonds(const Atom& atom, const Molecule& molecule, const Species& species,
                     AtomBondArray& bonds)
   {
      const int atomId  = int( &atom - &molecule.atom(0) );
      const Species::AtomBondIdArray bondIds = species.atomBondIds(atomId);
      const int nGroup = bondIds.size();
      const Bond* firstPtr = &molecule.bond(0);  // first Bond in molecule
      bonds.clear();
      for (int i = 0; i < nGroup; ++i) {
         bonds.append(firstPtr + bondIds[i]);
      }
   }

   /*
   * Fill an array of pointers to Bonds that contain an Atom.
   */
   void getAtomBonds(const Atom& atom, AtomBondArray& bonds)
   {
      const Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      getAtomBonds(atom, molecule, species, bonds);
   }
   #endif

} 
