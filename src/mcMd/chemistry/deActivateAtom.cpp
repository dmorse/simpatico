/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "deActivateAtom.h"
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
   * De-activate and atom an all associated groups.
   */
   void deActivateAtom(Atom& atom)
   {
      atom.setIsActive(false);

      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      int i;
      #ifdef INTER_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         const int nGroup = groupIds.size();
         for (i = 0; i < nGroup; ++i) {
            molecule.bond(groupIds[i]).setIsActive(false);
         }
      }
      #endif
      #ifdef INTER_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         const int nGroup = groupIds.size();
         for (i = 0; i < nGroup; ++i) {
            molecule.angle(groupIds[i]).setIsActive(false);
         }
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (species.nDihedral()) {
         const Species::AtomDihedralIdArray groupIds = species.atomDihedralIds(atomId);
         const int nGroup = groupIds.size();
         for (i = 0; i < nGroup; ++i) {
            molecule.dihedral(groupIds[i]).setIsActive(false);
         }
      }
      #endif
   }

   /*
   * Re-activate a temporarily de-activated atom and associated groups.
   *
   * Note: Each group is re-activated only if All atoms in the group are now active.
   */
   void reActivateAtom(Atom& atom)
   {
      atom.setIsActive(true);

      Molecule& molecule = atom.molecule();
      const Species& species = molecule.species();
      const int atomId = int( &atom - &molecule.atom(0) );
      int i, j;
      bool allActive;

      #ifdef INTER_BOND
      if (species.nBond()) {
         const Species::AtomBondIdArray groupIds = species.atomBondIds(atomId);
         const int nGroup = groupIds.size();
         Bond* groupPtr;
         for (i = 0; i < nGroup; ++i) {
            groupPtr = &molecule.bond(groupIds[i]);
            allActive = true;
            for (j = 0; j < 2; ++j) {
               if (!groupPtr->atom(j).isActive()) {
                  allActive = false;
                  break;
               }
            }
            if (allActive) {
               groupPtr->setIsActive(true);
            } else {
               groupPtr->setIsActive(false);
            }
         }
      }
      #endif

      #ifdef INTER_ANGLE
      if (species.nAngle()) {
         const Species::AtomAngleIdArray groupIds = species.atomAngleIds(atomId);
         const int nGroup = groupIds.size();
         Angle* groupPtr;
         for (i = 0; i < nGroup; ++i) {
            groupPtr = &molecule.angle(groupIds[i]);
            allActive = true;
            for (j = 0; j < 3; ++j) {
               if (!groupPtr->atom(j).isActive()) {
                  allActive = false;
                  break;
               }
            }
            if (allActive) {
               groupPtr->setIsActive(true);
            } else {
               groupPtr->setIsActive(false);
            }
         }
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (species.nDihedral()) {
         const Species::AtomDihedralIdArray groupIds = species.atomDihedralIds(atomId);
         const int nGroup = groupIds.size();
         Dihedral* groupPtr;
         for (i = 0; i < nGroup; ++i) {
            groupPtr = &molecule.dihedral(groupIds[i]);
            allActive = true;
            for (j = 0; j < 4; ++j) {
               if (!groupPtr->atom(j).isActive()) {
                  allActive = false;
                  break;
               }
            }
            if (allActive) {
               groupPtr->setIsActive(true);
            }
         }
      }
      #endif

   }

}
