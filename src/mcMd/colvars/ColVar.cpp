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
   * Output at end - empty default implementation.
   */
   void ColVar::unset()
   {  value_.unset(); }

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

}
