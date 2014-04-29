#ifndef UTIL_CONSTANTS_H
#define UTIL_CONSTANTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <complex>

namespace Util
{

   /**
   * Mathematical constants.
   *
   * \ingroup Math_Module
   */
   class Constants
   {

   public:

      /**
      * Initialize static constants.
      */
      static void initStatic();

      /**
      * Trigonometric constant Pi.
      */
      static const double               Pi;

      /**
      * Square root of -1.
      */
      static const std::complex<double> Im;

   };

} 
#endif
