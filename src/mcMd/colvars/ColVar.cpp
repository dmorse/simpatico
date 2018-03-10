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

   /*
   * Default constructor.
   */
   ColVar::ColVar() 
   {}

   /*
   * Default destructor.
   */
   ColVar::~ColVar()
   {}

   /*
   * Return value, compute if necesary.
   */
   double ColVar::value()
   {
      if (!value_.isSet()) {  
         compute();
      }
      return value_.value();
   }

   /*
   * Unset colvar value, i.e., mark as unknown.
   */
   void ColVar::unset()
   {  value_.unset(); }

   /*   
   * Return true iff value is set (i.e., known).
   */
   bool ColVar::isSet() const
   {  return value_.isSet(); }

}
