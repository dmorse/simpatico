#ifndef DDMD_SP_GROUP_STORAGE_H
#define DDMD_SP_GROUP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DSArray.h>       // member (template)
#include <ddMd/sp/chemistry/SpGroup.h>        // member (template argument)
#include <util/containers/ArrayIterator.h> // inline function

namespace DdMd 
{

   using namespace Util;

   /**
   * A container for SpGroup<N> objects.
   *
   * \ingroup DdMd_Sp_Storage_Module
   */
   template <int N>
   class SpGroupStorage : private DSArray< SpGroup<N> >
   {

   public:

      typedef ArrayIterator< SpGroup <N> > Iterator;

      typedef DSArray< SpGroup<N> > Base;
      using Base::allocate;
      using Base::begin;
      using Base::clear;
      using Base::size;
      using Base::capacity;

      /**
      * Destructor.
      */
      ~SpGroupStorage();

      /**
      * Append a new element to container and return its address.
      *
      * \return pointer to location of new SpGroup
      */
      SpGroup<N>* newPtr();

   private:

      using Base::resize;

   };

   // Function definition

   /*
   * Return pointer to location for new group, and add to container.
   */
   template <int N>
   SpGroup<N>* SpGroupStorage<N>::newPtr()
   {
      resize(size() + 1);
      return &(*this)[size()-1];
   }

   /*
   * Destructor
   */
   template <int N>
   SpGroupStorage<N>::~SpGroupStorage()
   {}

}
#endif
