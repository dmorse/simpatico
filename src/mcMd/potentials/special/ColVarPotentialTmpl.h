#ifndef MCMD_COLVAR_POTENTIAL_TMPL_H
#define MCMD_COLVAR_POTENTIAL_TMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/special/SpecialPotential.h> // base class

namespace McMd
{

   class System;

   using namespace Util;

   /**
   * SpecialPotential that is a function of a collective variable. 
   *
   * \ingroup McMd_Special_Module
   */
   template <class ColVarType, class BiasType>
   class ColVarPotentialTmpl : public SpecialPotential
   {

   public:

      /**
      * Constructor.
      */
      ColVarPotentialTmpl(System& system);

      /**
      * Destructor.
      */
      virtual ~ColVarPotentialTmpl();

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      /**
      * Compute total energy.
      */
      void computeEnergy();

      /**
      * Unset the total energy (mark as unknown or obsolete).
      */
      void unsetEnergy();

      /**
      * Add forces from this potential to all atomic forces.
      */
      void addForces();

   private:

      /// Pointer to parent System.
      System* systemPtr_;

      /// Collective Variable object
      ColVarType colVar_;

      /// Bias Function object
      BiasType bias_;

   };

}
#include "ColVarPotentialTmpl.tpp" 
#endif
