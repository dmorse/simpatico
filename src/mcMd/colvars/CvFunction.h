#ifndef MCMD_CV_FUNCTION_H
#define MCMD_CV_FUNCTION_H

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
   * Function of a collective variable (CV).
   *
   * This potential may be used to define a bias potential in MD or 
   * MC simulations.
   *
   * \ingroup McMd_CvFunction_Module
   */
   class CvFunction : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      CvFunction();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~CvFunction();

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
