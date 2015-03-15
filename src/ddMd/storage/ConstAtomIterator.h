#ifndef DDMD_CONST_ATOM_ITERATOR_H
#define DDMD_CONST_ATOM_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstPArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   class Atom;

   /**
   * Const iterator for all atoms owned by an AtomStorage.
   *
   * A ConstAtomIterator returns pointers const Atoms, and
   * thus does not allow atoms to be modified.
   *
   * \ingroup DdMd_Storage_Atom_Module
   */
   class ConstAtomIterator : public ConstPArrayIterator<Atom>
   {};

}
#endif
