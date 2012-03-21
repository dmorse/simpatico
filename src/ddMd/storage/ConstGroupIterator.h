#ifndef DDMD_CONST_GROUP_ITERATOR_H
#define DDMD_CONST_GROUP_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstPArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   template <int N> class Group; 

   /**
   * Iterator for all Group<N> objects owned by a System.
   */
   template <int N>
   class ConstGroupIterator : public ConstPArrayIterator< Group<N> >
   {};

   /**
   * Iterator for all Group<N> objects owned by a System.
   */
   template <int N>
   class ConstIncompleteGroupIterator : public ConstPArrayIterator< Group<N> >
   {};

}
#endif
