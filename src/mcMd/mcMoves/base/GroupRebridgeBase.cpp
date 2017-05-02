/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "GroupRebridgeBase.h"
#include <mcMd/mcSimulation/mc_potentials.h>
//#include <mcMd/potentials/pair/McPairPotential.h>
//#include <mcMd/potentials/bond/McBondPotential.h>
#include <mcMd/chemistry/getAtomGroups.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   GroupRebridgeBase::GroupRebridgeBase(McSystem& system) : 
      SystemMove(system)
   {
      // Precondition
      #ifdef SIMP_DIHEDRAL
      if (system.hasDihedralPotential()) {
         UTIL_THROW("GroupEndBase is unusable with dihedrals");
      }
      #endif
   } 
   
   /* 
   * Destructor.
   */
   GroupRebridgeBase::~GroupRebridgeBase() 
   {} 
  
   /*
   * Calculate the energy change of rebridging a tetrahedron formed from
   * 4 atoms.
   *
   *    a - b          a   b
   *             ->      X
   *    d - c          d   c
   *
   * Erase ab & cd bonds, and add ac & bd bonds. Energy change includes the
   * bonding potential changes and the pairing potential changes.
   */
   void GroupRebridgeBase::
   tetraEnergy(Atom* aPtr, Atom* bPtr, Atom* cPtr, Atom* dPtr,
               int bondType, double &energy)
   {
      Vector aPos = aPtr->position();
      Vector bPos = bPtr->position();
      Vector cPos = cPtr->position();
      Vector dPos = dPtr->position();
      double lab, lcd, lac, lbd;

      energy = 0.0;

      // Calculate bond lengthSq.
      lab = boundary().distanceSq(aPos, bPos);
      lcd = boundary().distanceSq(cPos, dPos);
      lac = boundary().distanceSq(aPos, cPos);
      lbd = boundary().distanceSq(bPos, dPos);

      // Subtract old bonding energies.
      energy -= system().bondPotential().energy(lab, bondType);
      energy -= system().bondPotential().energy(lcd, bondType);

      // Add new bonding energies.
      energy += system().bondPotential().energy(lac, bondType);
      energy += system().bondPotential().energy(lbd, bondType);
   
      #ifndef SIMP_NOPAIR
      int    aId = aPtr->typeId();
      int    bId = bPtr->typeId();
      int    cId = cPtr->typeId();
      int    dId = dPtr->typeId();

      // Subtract old pairing potentials.
      energy -= system().pairPotential().energy(lac, aId, cId);
      energy -= system().pairPotential().energy(lbd, bId, dId);

      // Add new pairing potentials.
      energy += system().pairPotential().energy(lab, aId, bId);
      energy += system().pairPotential().energy(lcd, cId, dId);
      #endif

      #ifdef SIMP_ANGLE
      if (system().hasAnglePotential()) {

         // Implementation assumed linear or ring molecular topology
         const Atom *amPtr(NULL), *bpPtr(NULL), *cmPtr(NULL), *dpPtr(NULL);
         AtomBondArray bonds;
         AtomAngleArray angles;
         const Bond *bondPtr;
         int   iBond, angleType;
   
         // find the neighboring atom pointers and angle type
   
         // angle type
         getAtomAngles(*aPtr, angles); 
         angleType = angles[0]->typeId();
   
         // amPtr
         getAtomBonds(*aPtr, bonds);
         for (iBond = 0; iBond < bonds.size(); ++iBond) {
            bondPtr = bonds[iBond];
            if (&bondPtr->atom(0) != aPtr && &bondPtr->atom(0) != bPtr) {
               amPtr = &bondPtr->atom(0);
            } else if (&bondPtr->atom(1) != aPtr && &bondPtr->atom(1) != bPtr){
               amPtr = &bondPtr->atom(1);
            }
         }
   
         // bpPtr
         getAtomBonds(*bPtr, bonds);
         for (iBond = 0; iBond < bonds.size(); ++iBond) {
            bondPtr = bonds[iBond];
            if (&bondPtr->atom(0) != bPtr && &bondPtr->atom(0) != aPtr) {
               bpPtr = &bondPtr->atom(0);
            } else if (&bondPtr->atom(1) != bPtr && &bondPtr->atom(1) != aPtr){
               bpPtr = &bondPtr->atom(1);
            }
         }
   
         // cmPtr
         getAtomBonds(*cPtr, bonds);
         for (iBond = 0; iBond < bonds.size(); ++iBond) {
            bondPtr = bonds[iBond];
            if (&bondPtr->atom(0) != cPtr && &bondPtr->atom(0) != dPtr) {
               cmPtr = &bondPtr->atom(0);
            } else if (&bondPtr->atom(1) != cPtr && &bondPtr->atom(1) != dPtr){
               cmPtr = &bondPtr->atom(1);
            }
         }
   
         // dpPtr
         getAtomBonds(*dPtr, bonds);
         for (iBond = 0; iBond < bonds.size(); ++iBond) {
            bondPtr = bonds[iBond];
            if (&bondPtr->atom(0) != dPtr && &bondPtr->atom(0) != cPtr) {
               dpPtr = &bondPtr->atom(0);
            } else if (&bondPtr->atom(1) != cPtr && &bondPtr->atom(1) != dPtr){
               dpPtr = &bondPtr->atom(1);
            }
         }
   
         // subtract energies of angles involving bond a-b
         energy -= angleEnergy(*amPtr, *aPtr, *bPtr, angleType);
         energy -= angleEnergy(*aPtr, *bPtr, *bpPtr, angleType);
   
         // subtract energies of angles involving bond c-d
         energy -= angleEnergy(*cmPtr, *cPtr, *dPtr, angleType);
         energy -= angleEnergy(*cPtr, *dPtr, *dpPtr, angleType);
   
         // add energies of angles involving bond a-c
         energy += angleEnergy(*amPtr, *aPtr, *cPtr, angleType);
         energy += angleEnergy(*aPtr, *cPtr, *cmPtr, angleType);
   
         // add energies of angles involving bond b-d
         energy += angleEnergy(*bpPtr, *bPtr, *dPtr, angleType);
         energy += angleEnergy(*bPtr, *dPtr, *dpPtr, angleType);
   
      }
      #endif
   }

   /*
   * Calculate the energy change of rebridging an octahedron formed
   * from 6 atoms.
   *
   *    a - m - b          a   m   b
   *                 ->      X   X
   *    c - n - d          c   n   d
   *
   * Erase 4 bonds: am, mb, cn, nd
   * Add   4 bonds: an, nb, cm, md
   * The erased bonds contribute to the pairing interactions, but
   * the pair of atoms involved in the newly added bonds lost their pairing
   * interaction.
   */
   void GroupRebridgeBase::
   octaEnergy(Atom* aPtr, Atom* bPtr, Atom* cPtr, Atom* dPtr,
              Atom* mPtr, Atom* nPtr, int bondType, double &energy)
   {
      Vector aPos = aPtr->position();
      Vector bPos = bPtr->position();
      Vector cPos = cPtr->position();
      Vector dPos = dPtr->position();
      Vector mPos = mPtr->position();
      Vector nPos = nPtr->position();
      double lam, lmb, lcn, lnd;
      double lan, lnb, lcm, lmd;

      energy = 0.0;

      // Calculate atom pair distances squared.
      lam = boundary().distanceSq(aPos, mPos);
      lmb = boundary().distanceSq(mPos, bPos);
      lcn = boundary().distanceSq(cPos, nPos);
      lnd = boundary().distanceSq(nPos, dPos);

      lan = boundary().distanceSq(aPos, nPos);
      lnb = boundary().distanceSq(nPos, bPos);
      lcm = boundary().distanceSq(cPos, mPos);
      lmd = boundary().distanceSq(mPos, dPos);

      // Subtract old bonding energies.
      energy -= system().bondPotential().energy(lam, bondType);
      energy -= system().bondPotential().energy(lmb, bondType);
      energy -= system().bondPotential().energy(lcn, bondType);
      energy -= system().bondPotential().energy(lnd, bondType);

      // Add new bonding energies.
      energy += system().bondPotential().energy(lan, bondType);
      energy += system().bondPotential().energy(lnb, bondType);
      energy += system().bondPotential().energy(lcm, bondType);
      energy += system().bondPotential().energy(lmd, bondType);
  
      #ifndef SIMP_NOPAIR
      int    aId = aPtr->typeId();
      int    bId = bPtr->typeId();
      int    cId = cPtr->typeId();
      int    dId = dPtr->typeId();
      int    mId = mPtr->typeId();
      int    nId = nPtr->typeId();

      // Subtract old pairing potentials.
      energy -= system().pairPotential().energy(lan, aId, nId);
      energy -= system().pairPotential().energy(lnb, nId, bId);
      energy -= system().pairPotential().energy(lcm, cId, mId);
      energy -= system().pairPotential().energy(lmd, mId, dId);

      // Add new pairing potentials.
      energy += system().pairPotential().energy(lam, aId, mId);
      energy += system().pairPotential().energy(lmb, mId, bId);
      energy += system().pairPotential().energy(lcn, cId, nId);
      energy += system().pairPotential().energy(lnd, nId, dId);
      #endif
   }
}
