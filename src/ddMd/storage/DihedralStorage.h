#ifndef DDMD_DIHEDRAL_STORAGE_H
#define DDMD_DIHEDRAL_STORAGE_H

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
   * Container for Group<4> (dihedral) objects.
   *
   * Only differences from GroupStorage<4> is the className,
   * which is set in the constructor and effects the param 
   * file format.
   *
   * \ingroup DdMd_Storage_Module
   */
   class DihedralStorage : public GroupStorage<4>
   {

   public:

      /**
      * Constructor
      */
      DihedralStorage();

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      *   
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

   };

}
#endif
