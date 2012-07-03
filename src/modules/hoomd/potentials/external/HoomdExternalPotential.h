#ifdef INTER_EXTERNAL
#ifndef HOOMD_EXTERNAL_POTENTIAL_H
#define HOOMD_EXTERNAL_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   /**
   * Abstract interface common to all HOOMD external potentials
   * (used for querying the evaluator name)
   */
   class HoomdExternalPotential
   {

   public:
      
      /**
      * get the internal name of the HOOMD evaluator
      */ 
      virtual std::string hoomdName() const = 0;

   };

}

#endif
#endif
