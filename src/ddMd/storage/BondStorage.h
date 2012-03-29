#ifndef DDMD_BOND_STORAGE_H
#define DDMD_BOND_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorage.h"

namespace DdMd
{

   /**
   * AtomStorage for bonds.
   *
   * \ingroup DdMd_Storage_Module
   */
   class BondStorage : public GroupStorage<2>
   {

   public:

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      */
      virtual void readParam(std::istream& in);

   };

}
#endif
