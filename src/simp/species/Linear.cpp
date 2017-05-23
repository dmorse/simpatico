/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Linear.h"

#include <sstream>

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   Linear::Linear() 
   : Species()
     #ifdef SIMP_ANGLE
     , hasAngles_(0)
     #endif
     #ifdef SIMP_DIHEDRAL
     , hasDihedrals_(0)
     #endif
   {}

   /*
   * Destructor.
   */
   Linear::~Linear() 
   {}

   /*
   * Build the chemical structure of a flexible chain molecule.
   */
   void Linear::buildLinear() 
   {

      // Preconditions
      if (nAtom() < 2) {
         UTIL_THROW("nAtom < 2");
      }
      if (nBond() != nAtom() - 1) {
         UTIL_THROW("nBond != nAtom - 1");
      }
      #ifdef SIMP_ANGLE
      if (hasAngles_) {
         if (nAngle() != nAtom() - 2) {
            UTIL_THROW("nAngle != nAtom - 2");
         }
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals_) {
         if (nAtom() > 3) {
            if (nDihedral() != nAtom() - 3) {
               UTIL_THROW("nDihedral != nAtom - 3");
            }
         } else {
            if (nDihedral() != 0) {
               UTIL_THROW("nDihedral != 0 when nAtom < 4");
            }
         }
      }
      #endif

      int i; 

      allocate();

      // Set Atom Type Ids
      for (i = 0; i < nAtom(); ++i) {
         setAtomType(i, calculateAtomTypeId(i));
      }

      // Build Bonds
      for (i = 0; i < nBond(); ++i) {
         makeBond(i, i, i+1, calculateBondTypeId(i));
      }

      #ifdef SIMP_ANGLE
      // Build Angles 
      if (hasAngles_) {
         for (i = 0; i < nAngle(); ++i) {
            makeAngle(i, i, i+1, i+2, calculateAngleTypeId(i));
         }
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Build Dihedrals.
      if (hasDihedrals_) {
         for (i = 0; i < nDihedral(); ++i) {
            makeDihedral(i, i, i+1, i+2, i+3, calculateDihedralTypeId(i));
         }
      }
      #endif

   }
   
}
