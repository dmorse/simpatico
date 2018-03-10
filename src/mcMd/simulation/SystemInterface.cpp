/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// namespace McMd
#include "SystemInterface.h"
#include "Simulation.h"
#include <simp/species/Species.h>

#include <fstream>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   SystemInterface::SystemInterface(System& parent)
    : simulationPtr_(parent.simulationPtr_),
      systemPtr_(&parent),
      moleculeSetsPtr_(parent.moleculeSetsPtr_),
      boundaryPtr_(parent.boundaryPtr_)
      #ifdef SIMP_BOND
      , hasBonds_(parent.simulation().nBondType() > 0)
      #endif
      #ifdef SIMP_ANGLE
      , hasAngles_(parent.simulation().nAngleType() > 0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , hasDihedrals_(parent.simulation().nDihedralType() > 0)
      #endif
      #ifdef MCMD_LINK
      , hasLinks_(parent.simulation().nLinkType() > 0)
      #endif
      #ifdef SIMP_EXTERNAL
      , hasExternal_(parent.simulation().hasExternal())
      #endif
      #ifdef SIMP_TETHER
      , hasTethers_(parent.simulation().hasTether())
      #endif
   {}

   /*
   * Destructor.
   */
   SystemInterface::~SystemInterface()
   {}

   /*
   * Is this System empty (i.e., devoid of Molecules) ?
   */
   bool SystemInterface::isEmpty() const
   {
      for (int i = 0; i < simulation().nSpecies(); ++i) {
         if (nMolecule(i) != 0) return false;
      }
      return true;
   }

   /*
   * Return the total number of atoms in this System.
   */
   int SystemInterface::nAtom() const
   {
      int sum = 0;
      for (int i = 0; i < simulation().nSpecies(); ++i) {
         sum += nMolecule(i)*simulation().species(i).nAtom();
      }
      return sum;
   }

}
