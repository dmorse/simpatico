/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Colvar.h"

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   Colvar::Colvar() 
   {}

   /*
   * Default destructor.
   */
   Colvar::~Colvar()
   {}

   /*
   * Output at end - empty default implementation.
   */
   void Colvar::unset()
   {  value_.unset(); }

   /*
   * Return value, compute if necesary.
   */
   double Colvar::value()
   {
      if (!value_.isSet()) {  
         compute();
      }
      return value_.value();
   }

}
