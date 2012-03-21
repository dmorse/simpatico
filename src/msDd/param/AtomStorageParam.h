#ifndef MSDD_ATOM_STORAGE_PARAM_H
#define MSDD_ATOM_STORAGE_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class

namespace MsDd 
{

   using namespace Util;

   /**
   * Parameters for an DdMd::AtomStorageParam.
   */
   class AtomStorageParam : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      AtomStorageParam();

      /**
      * Destructor.
      */
      ~AtomStorageParam();

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - atomCapacity      [int]  max number of atoms owned by processor.
      *  - ghostCapacity     [int]  max number of ghosts on this processor.
      *  - totalatomCapacity [int]  max number of atoms on all processors.
      *
      * \param in input parameter stream.
      */
      virtual void readParam(std::istream& in);

      /**
      * Return capacity for atoms (excluding ghosts).
      */
      int atomCapacity() const;

      /**
      * Return capacity for ghost atoms
      */
      int ghostCapacity() const;

      /**
      * Return maximum number of atoms on all processors.
      *
      * Atom ids are labelled from 0, ..., totalAtomCapacity-1
      */
      int totalAtomCapacity() const;


   private:

      // Capacity for local atoms on this processor.
      int atomCapacity_;

      // Capacity for ghost atoms on this processors.
      int ghostCapacity_;

      // Maximum number of atoms on all processors, maximum id + 1
      int totalAtomCapacity_;

   };

   // Inline method definitions

   inline int AtomStorageParam::atomCapacity() const
   { return atomCapacity_; }

   inline int AtomStorageParam::ghostCapacity() const
   { return ghostCapacity_; }

   inline int AtomStorageParam::totalAtomCapacity() const
   { return totalAtomCapacity_; }

}
#endif
