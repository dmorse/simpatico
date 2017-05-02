/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "getAtomGroups.h"
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Molecule.h>

namespace McMd
{

   using namespace Util;

   #ifdef SIMP_BOND
   /*
   * Fill an array of pointers to Bonds that contain an Atom.
   */
   void getAtomBonds(const Atom& atom, AtomBondArray& groups)
   {
      groups.clear();
      const Molecule& molecule = atom.molecule();
      if (molecule.nBond()) {
         const Species& species = molecule.species();
         const int atomId  = int( &atom - &molecule.atom(0) );
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         const int nGroup = groupIds.size();
         const Bond* firstPtr = &molecule.bond(0);  // first group in molecule
         for (int i = 0; i < nGroup; ++i) {
            groups.append(firstPtr + groupIds[i]);
         }
      }
   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Fill an array of pointers to Angle objects that contain an atom.
   */
   void getAtomAngles(const Atom& atom, AtomAngleArray& groups)
   {
      groups.clear();
      const Molecule& molecule = atom.molecule();
      if (molecule.nAngle()) {
         const Species& species = molecule.species();
         const int atomId  = int( &atom - &molecule.atom(0) );
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         const int nGroup = groupIds.size();
         const Angle* firstPtr = &molecule.angle(0);  // first group in molecule
         for (int i = 0; i < nGroup; ++i) {
            groups.append(firstPtr + groupIds[i]);
         }
      }
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Fill an array of pointers to Dihedrals that contain an Atom.
   */
   void getAtomDihedrals(const Atom& atom, AtomDihedralArray& groups)
   {
      groups.clear();
      const Molecule& molecule = atom.molecule();
      if (molecule.nDihedral()) {
         const Species& species = molecule.species();
         const int atomId  = int( &atom - &molecule.atom(0) );
         const Species::AtomDihedralIdArray groupIds = species.atomDihedralIds(atomId);
         const int nGroup = groupIds.size();
         const Dihedral* firstPtr = &molecule.dihedral(0);  // first group in molecule
         for (int i = 0; i < nGroup; ++i) {
            groups.append(firstPtr + groupIds[i]);
         }
      }
   }
   #endif

}
