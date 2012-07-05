#ifndef FEQ_H
#define FEQ_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <cmath>

namespace Util
{

   /**
   * Are two floating point numbers equal to within round-off error?
   * 
   * \ingroup Math_Module
   */
   inline bool feq(double x, double y, double eps = 1.0E-12)
   {  return (fabs(x-y) < eps); }

}
#endif
