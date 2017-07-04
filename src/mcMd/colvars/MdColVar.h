#ifndef MCMD_MD_COLVAR_H
#define MCMD_MD_COLVAR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ColVar.h"

namespace McMd
{

   using namespace Util;

   /**
   * Collective variables for MD simulation, with atom forces.
   *
   * Added addForces() function computes forces on atoms.
   *
   * \ingroup McMd_ColVar_Module
   */
   class MdColVar : public ColVar
   {

   public:

      /**
      * Constructor.
      */
      MdColVar();

      /**
      * Destructor.
      */
      virtual ~MdColVar();

      /**
      * Add atom forces arising from bias potential.
      * 
      * The force on each atom is computed as -1 times the product 
      * of the derivative d(Cv)/dr of the collective variable Cv 
      * with respect to atom position r and the derivative dW/d(Cv) 
      * of the 1D bias potential W(Cv) with respect to the colvar
      * Cv, i.e., f = - d(Cv)/dr * dW/d(Cv).
      *
      * \param dWdCv derivative of bias potential W w/r to colvar.
      */
      virtual void addForces(double dWdCv) = 0;

   };

}
#endif
