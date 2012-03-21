#ifndef MSDD_GROUP_STORAGE_PARAM_H
#define MSDD_GROUP_STORAGE_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/global.h>


namespace MsDd
{

   using namespace Util;

   /**
   * Parameters for a DdMd::GroupStorage<N> container.
   */
   template <int N>
   class GroupStorageParam : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      GroupStorageParam();

      /**
      * Destructor.
      */
      ~GroupStorageParam();

      /**
      * Read parameters, allocate memory and initialize.
      *
      * Parameters (file format):
      *  - capacity      [int]  max number of groups owned by processor.
      *  - totalcapacity [int]  max number of groups on all processors.
      */
      virtual void readParam(std::istream& in);

      int capacity() const;

      int totalCapacity() const;

   private:

      // Capacity for local groups on this processor.
      int capacity_;

      // Maximum number of groups on all processors, maximum id + 1
      int totalCapacity_;

   };

   // Inline method definitions

   template <int N>
   inline int GroupStorageParam<N>::capacity() const
   {  return capacity_; }

   template <int N>
   inline int GroupStorageParam<N>::totalCapacity() const
   {  return totalCapacity_; }

   // Non-inline method templates.

   /*
   * Default constructor.
   */
   template <int N>
   GroupStorageParam<N>::GroupStorageParam()
    : capacity_(0),
      totalCapacity_(0)
   {}
 
   /*
   * Destructor.
   */
   template <int N>
   GroupStorageParam<N>::~GroupStorageParam()
   {}

   /*
   * Read parameters.
   */
   template <int N>
   void GroupStorageParam<N>::readParam(std::istream& in)
   {
      read<int>(in, "capacity", capacity_);
      read<int>(in, "totalCapacity", totalCapacity_);
   }

}
#endif
