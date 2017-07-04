#ifndef MCMD_CV_POTENTIAL_H
#define MCMD_CV_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>

namespace McMd
{

   using namespace Util;

   /**
   * Bias potential for a collective variable (CV).
   *
   * This potential may be used in MD as well as MC simulations.
   *
   * \ingroup McMd_CvPotential_Module
   */
   class CvPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      CvPotential();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~CvPotential();

      /**
      * Compute and return the bias potential value.
      *
      * \param colVar collective variable value. 
      */
      virtual double value(double colVar) = 0;

      /**
      * Compute and return the derivative of the bias potential.
      *
      * \param colVar collective variable value. 
      */
      virtual double derivative(double colVar) = 0;

   };

}
#endif
