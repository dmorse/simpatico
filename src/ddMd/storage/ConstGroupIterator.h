#ifndef DDMD_CONST_GROUP_ITERATOR_H
#define DDMD_CONST_GROUP_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstPArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   template <int N> class Group; 

   /**
   * Const iterator for all Group < N > objects owned by a GroupStorage < N >.
   *
   * This iterator provides read-only (const) access to Group < N > objects.
   *
   * \ingroup DdMd_Storage_Module
   */
   template <int N>
   class ConstGroupIterator : public ConstPArrayIterator< Group<N> >
   {};

}
#endif
