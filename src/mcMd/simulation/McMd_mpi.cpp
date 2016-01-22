/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMd_mpi.h"
#ifdef UTIL_MPI

namespace Util
{
   /*
   * Initialize MpiTraits< SpeciesGroup<4> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<4> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<4> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<3> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<3> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<3> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<2> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<2> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<2> >::hasType = false;

   /*
   * Initialize MpiTraits< Pair<int> >
   */
   MPI::Datatype MpiTraits< Pair<int> >::type = MPI::BYTE;
   bool MpiTraits< Pair<int> >::hasType = false;

}

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <mcMd/species/SpeciesGroup.tpp>
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
      if (!Util::MpiTraits< SpeciesGroup<2> >::hasType) {
         SpeciesGroup<2>::commitMpiType();
      }
      if (!Util::MpiTraits< SpeciesGroup<3> >::hasType) {
         SpeciesGroup<3>::commitMpiType();
      }
      if (!Util::MpiTraits< SpeciesGroup<4> >::hasType) {
         SpeciesGroup<4>::commitMpiType();
      }
      if (!Util::MpiTraits< PairSelector >::hasType) {
         PairSelector::commitMpiType();
      }
   }

}
#endif // ifdef  UTIL_MPI
