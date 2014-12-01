/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeActivator.h"
#include <mcMd/species/Species.h>     
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Molecule.h>
#ifdef INTER_BOND
#include <mcMd/chemistry/Bond.h>    
#endif
#ifdef INTER_ANGLE
#include <mcMd/chemistry/Angle.h>   
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/chemistry/Dihedral.h>  
#endif

namespace McMd
{

   using namespace Util;

   /*
   * De-activate an atom and all associated groups.
   */
   void DeActivator::deActivate(Atom& atom)
   {
      assert(atom.isActive());
      atom.setIsActive(false);

      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      int i;
      #ifdef INTER_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         for (i = 0; i < groupIds.size(); ++i) {
            molecule.bond(groupIds[i]).incrementInactive();
            assert(molecule.bond(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef INTER_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         for (i = 0; i < groupIds.size(); ++i) {
            molecule.angle(groupIds[i]).incrementInactive();
            assert(molecule.angle(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (species.nDihedral()) {
         const Species::AtomDihedralIdArray groupIds = species.atomDihedralIds(atomId);
         for (i = 0; i < groupIds.size(); ++i) {
            molecule.dihedral(groupIds[i]).incrementInactive();
            assert(molecule.dihedral(groupIds[i]).checkInactive());
         }
      }
      #endif
   }

   /*
   * Re-activate a temporarily de-activated atom and associated groups.
   *
   * Note: Each group is re-activated only if All atoms in the group are now active.
   */
   void DeActivator::reActivate(Atom& atom)
   {
      assert(!atom.isActive());
      atom.setIsActive(true);

      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      assert(atomId >= 0);
      assert(atomId < molecule.nAtom());

      #ifdef INTER_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         for (int i = 0; i < groupIds.size(); ++i) {
            molecule.bond(groupIds[i]).decrementInactive();
            assert(molecule.bond(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef INTER_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         for (int i = 0; i < groupIds.size(); ++i) {
            molecule.angle(groupIds[i]).decrementInactive();
            assert(molecule.angle(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (species.nDihedral()) {
         const Species::AtomDihedralIdArray groupIds = species.atomDihedralIds(atomId);
         for (int i = 0; i < groupIds.size(); ++i) {
            molecule.dihedral(groupIds[i]).decrementInactive();
            assert(molecule.dihedral(groupIds[i]).checkInactive());
         }
      }
      #endif

   }

   /*
   * Activate all atoms and groups in this molecule.
   */
   void DeActivator::activate(Molecule& molecule)
   {
      int i;
      for (i = 0; i < molecule.nAtom(); ++i) {
         molecule.atom(i).setIsActive(true);
      }
      #ifdef INTER_BOND
      for (i = 0; i < molecule.nBond(); ++i) {
         molecule.bond(i).activate();
      }
      #endif
      #ifdef INTER_ANGLE
      for (i = 0; i < molecule.nAngle(); ++i) {
         molecule.angle(i).activate();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      for (i = 0; i < molecule.nDihedral(); ++i) {
        molecule.dihedral(i).activate();
      }
      #endif
   }

}
