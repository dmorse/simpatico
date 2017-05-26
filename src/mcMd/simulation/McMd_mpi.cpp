/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#ifdef UTIL_MPI
#include "McMd_mpi.h"
#include <simp/species/SpeciesGroup.h>

namespace Util
{

   /*
   * Initialize MpiTraits< Pair<int> >
   */
   MPI::Datatype MpiTraits< Pair<int> >::type = MPI::BYTE;
   bool MpiTraits< Pair<int> >::hasType = false;

}

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <simp/species/SpeciesGroup.tpp>
#include <mcMd/analyzers/util/PairSelector.h>

namespace McMd
{

   void commitMpiTypes()
   {
      if (!Util::MpiTraits<Util::Vector>::hasType) {
         Util::Vector::commitMpiType();
      }
      if (!Util::MpiTraits<Util::IntVector>::hasType) {
         Util::IntVector::commitMpiType();
      }
      if (!Util::MpiTraits< Util::Pair<int> >::hasType) {
         Util::Pair<int>::commitMpiType();
      }
      if (!Util::MpiTraits< PairSelector >::hasType) {
         PairSelector::commitMpiType();
      }
      Simp::commitMpiSpeciesGroupTypes();
   }

}
#endif // ifdef  UTIL_MPI
