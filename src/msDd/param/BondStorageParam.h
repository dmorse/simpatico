#ifndef MSDD_BOND_STORAGE_PARAM_H
#define MSDD_BOND_STORAGE_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorageParam.h"

namespace MsDd 
{

   using namespace Util;

   /**
   * Parameters for an DdMd::BondStorageParam.
   */
   class BondStorageParam : public GroupStorageParam<2>
   {

   public:

      /**
      * Constructor.
      */
      BondStorageParam();

      /**
      * Destructor.
      */
      ~BondStorageParam();

      /**
      * Read parameters.
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

   };

}
#endif
