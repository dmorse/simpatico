#ifndef DDMD_ATOM_ITERATOR_H
#define DDMD_ATOM_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArrayIterator.h>

namespace DdMd
{

   using namespace Util;
   class Atom;

   /**
   * Iterator for all atoms owned by an AtomStorage.
   *
   * \ingroup DdMd_Storage_Module
   */
   class AtomIterator : public PArrayIterator<Atom>
   {};

}
#endif
