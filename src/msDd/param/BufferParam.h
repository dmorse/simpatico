#ifndef MSDD_BUFFER_PARAM_H
#define MSDD_BUFFER_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/
#include <util/param/ParamComposite.h>  // base class

namespace MsDd
{

   using namespace Util;

   /**
   * BufferParam for sending local or ghost atoms.
   */
   class BufferParam: public ParamComposite 
   {

   public:

      /**
      * Constructor.
      */
      BufferParam();

      /**
      * Destructor (empty).
      */
      virtual ~BufferParam();

      /**
      * Read capacities.
      */
      void readParam(std::istream& in);

      /**
      * Maximum number of atoms for which space is available.
      */
      int atomCapacity() const;

      /**
      * Maximum number of ghost atoms for which space is available.
      */
      int ghostCapacity() const;

   private:

      /// Maximum number of local atoms in buffer.
      int atomCapacity_;

      /// Maximum number of ghost atoms in buffer.
      int ghostCapacity_;

   };

}
#endif
