#ifndef SPAN_GROUP_STORAGE_H
#define SPAN_GROUP_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DSArray.h>       // member (template)
#include <spAn/chemistry/Group.h>        // member (template argument)
#include <util/containers/ArrayIterator.h> // inline function

namespace SpAn 
{

   using namespace Util;

   /**
   * A container for Group<N> objects.
   *
   * \ingroup SpAn_Storage_Module
   */
   template <int N>
   class GroupStorage : private DSArray< Group<N> >
   {

   public:

      typedef ArrayIterator< Group <N> > Iterator;

      typedef DSArray< Group<N> > Base;
      using Base::allocate;
      using Base::begin;
      using Base::clear;
      using Base::size;
      using Base::capacity;

      /**
      * Destructor.
      */
      ~GroupStorage();

      /**
      * Append a new element to container and return its address.
      *
      * \return pointer to location of new Group
      */
      Group<N>* newPtr();

   private:

      using Base::resize;

   };

   // Function definition

   /*
   * Return pointer to location for new group, and add to container.
   */
   template <int N>
   Group<N>* GroupStorage<N>::newPtr()
   {
      resize(size() + 1);
      return &(*this)[size()-1];
   }

   /*
   * Destructor
   */
   template <int N>
   GroupStorage<N>::~GroupStorage()
   {}

}
#endif
