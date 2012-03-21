#ifdef  INTER_TETHER
#ifndef MCMD_MD_TETHER_POTENTIAL_H
#define MCMD_MD_TETHER_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/tether/TetherPotential.h> // base class

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * A TetherPotential for MD simulation (abstract).
   *
   * \ingroup Tether_Module
   */
   class MdTetherPotential : public TetherPotential 
   {

   #if 0
   public:

      /**   
      * Constructor.
      */
      MdTetherPotential();

      /** 
      * Destructor.
      */
      virtual ~MdTetherPotential();

      /**
      */
      virtual void addForces() = 0;
   #endif

   };

} 
#endif
#endif
