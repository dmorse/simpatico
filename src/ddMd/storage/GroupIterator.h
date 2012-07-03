#ifndef DDMD_GROUP_ITERATOR_H
#define DDMD_GROUP_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   template <int N> class Group; 

   /**
   * Iterator for all Group < N > objects owned by a GroupStorage< N >.
   *
   * \ingroup DdMd_Storage_Module
   */
   template <int N>
   class GroupIterator : public PArrayIterator< Group<N> >
   {};

}
#endif
