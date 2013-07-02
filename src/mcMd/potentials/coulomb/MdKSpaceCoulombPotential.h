#ifndef MCMD_MD_KSPACE_COULOMB_POTENTIAL_H
#define MCMD_MD_KSPACE_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdKSpaceCoulombPotential.h"

namespace McMd
{

   using namespace Util;

   /**
   * Long-range part of Coulomb potential for MD.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdKSpaceCoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor .
      */
      MdKSpaceCoulombPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~MdKSpaceCoulombPotential();

      /**
      * Add bond forces to all atomic forces.
      */
      void addForces();

   };

} 
#endif
