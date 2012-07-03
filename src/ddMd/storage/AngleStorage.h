#ifndef DDMD_ANGLE_STORAGE_H
#define DDMD_ANGLE_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "GroupStorage.h"

namespace DdMd
{

   /**
   * Container for angles objects.
   *
   * \ingroup DdMd_Storage_Module
   */
   class AngleStorage : public GroupStorage<3>
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
