#ifndef MCMD_MD_COULOMB_POTENTIAL_H
#define MCMD_MD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdCoulombPotential.h"

namespace McMd
{

   using namespace Util;

   /**
   * Long-range part of Coulomb potential for MD.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdCoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor .
      */
      MdCoulombPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~MdCoulombPotential();

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      void addKSpaceForces();

   };

} 
#endif
