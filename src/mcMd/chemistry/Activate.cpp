/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Activate.h"
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Molecule.h>
#include <simp/species/Species.h>     
#ifdef SIMP_BOND
#include <mcMd/chemistry/Bond.h>    
#endif
#ifdef SIMP_ANGLE
#include <mcMd/chemistry/Angle.h>   
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/chemistry/Dihedral.h>  
#endif

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * De-activate an atom and update associated groups.
   */
   void Activate::deactivate(Atom& atom)
   {
      // Precondition: Atom must be active
      if (!atom.isActive()) {
         UTIL_THROW("Atom is not active");
      }

      // De-activate atom
      atom.setIsActive(false);

      // Update associated groups
      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      int i;
      #ifdef SIMP_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         for (i = 0; i < groupIds.size(); ++i) {
            molecule.bond(groupIds[i]).incrementInactive();
            assert(molecule.bond(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         for (i = 0; i < groupIds.size(); ++i) {
            molecule.angle(groupIds[i]).incrementInactive();
            assert(molecule.angle(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
   * Re-activate a temporarily de-activated atom and update associated groups.
   */
   void Activate::reactivate(Atom& atom)
   {
      // Precondition: Atom must be inactive
      if (atom.isActive()) {
         UTIL_THROW("Atom already active");
      }

      // Re-activate atom
      atom.setIsActive(true);

      // Update associated groups
      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      assert(atomId >= 0);
      assert(atomId < molecule.nAtom());
      #ifdef SIMP_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         for (int i = 0; i < groupIds.size(); ++i) {
            molecule.bond(groupIds[i]).decrementInactive();
            assert(molecule.bond(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         for (int i = 0; i < groupIds.size(); ++i) {
            molecule.angle(groupIds[i]).decrementInactive();
            assert(molecule.angle(groupIds[i]).checkInactive());
         }
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
   void Activate::activate(Molecule& molecule)
   {
      int i;
      for (i = 0; i < molecule.nAtom(); ++i) {
         molecule.atom(i).setIsActive(true);
      }
      #ifdef SIMP_BOND
      for (i = 0; i < molecule.nBond(); ++i) {
         molecule.bond(i).activate();
      }
      #endif
      #ifdef SIMP_ANGLE
      for (i = 0; i < molecule.nAngle(); ++i) {
         molecule.angle(i).activate();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      for (i = 0; i < molecule.nDihedral(); ++i) {
        molecule.dihedral(i).activate();
      }
      #endif
   }

}
