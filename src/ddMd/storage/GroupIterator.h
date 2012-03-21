#ifndef DDMD_GROUP_ITERATOR_H
#define DDMD_GROUP_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   template <int N> class Group; 

   /**
   * Iterator for all Group<N> objects owned by this processor.
   */
   template <int N>
   class GroupIterator : public PArrayIterator< Group<N> >
   {};

   /**
   * Iterator for all incomplete Group<N> objects owned by this processor.
   */
   template <int N>
   class IncompleteGroupIterator : public PArrayIterator< Group<N> >
   {};

}
#endif
