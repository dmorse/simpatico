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
   * Container for Group<3> (angle) objects.
   *
   * Only differences from GroupStorage<3> is the className,
   * which is set in the constructor and effects the param 
   * file format.
   *
   * \ingroup DdMd_Storage_Module
   */
   class AngleStorage : public GroupStorage<3>
   {

   public:

      /**
      * Constructor
      */
      AngleStorage();

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      */
      virtual void readParameters(std::istream& in);

   };

}
#endif
