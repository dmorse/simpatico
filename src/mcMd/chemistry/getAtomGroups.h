#ifndef MCMD_GET_ATOM_GROUPS_H
#define MCMD_GET_ATOM_GROUPS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/FSArray.h>      // used in typedefs
#ifdef SIMP_BOND
#include <mcMd/chemistry/Bond.h>          // typedef
#endif
#ifdef SIMP_ANGLE
#include <mcMd/chemistry/Angle.h>         // typedef
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/chemistry/Dihedral.h>      // typedef
#endif
#include <simp/species/Species.h>         // used in typedefs

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   class Atom;
   class Molecule;

   #ifdef SIMP_BOND
   /**
   * Array to hold pointers to bonds that contain a specific atom.
   */
   typedef FSArray<const Bond*, Species::MaxBondPerAtom> 
           AtomBondArray;

   /**
   * Fill an array of pointers to Bonds that contain an Atom.
   *
   * \param atom the Atom of interest
   * \param bonds an array to fill with pointers to bonds
   */
   void getAtomBonds(const Atom& atom, AtomBondArray& bonds);
   #endif

   #ifdef SIMP_ANGLE
   /**
   * Array to hold pointers to angles that contain a specific atom.
   */
   typedef FSArray<const Angle*, Simp::Species::MaxAnglePerAtom> 
           AtomAngleArray;

   /**
   * Fill an array of pointers to Angles that contain an Atom.
   *
   * \param atom the Atom of interest
   * \param angles an array to fill with pointers to angles
   */
   void getAtomAngles(const Atom& atom, AtomAngleArray& angles);
   #endif

   #ifdef SIMP_DIHEDRAL
   /**
   * Array to hold pointers to Dihedrals that contain a specific atom.
   */
   typedef FSArray<const Dihedral*, Simp::Species::MaxDihedralPerAtom> 
           AtomDihedralArray;

   /**
   * Fill an array of pointers to Dihedrals that contain an Atom.
   *
   * \param atom the Atom of interest
   * \param dihedrals an array to fill with pointers to dihedrals
   */
   void 
   getAtomDihedrals(const Atom& atom, AtomDihedralArray& dihedrals);
   #endif

} 
#endif
