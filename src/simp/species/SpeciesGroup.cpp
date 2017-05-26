#ifdef UTIL_MPI

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesGroup.h"
#include "SpeciesGroup.tpp"

namespace Util
{
   /*
   * Initialize MpiTraits< SpeciesGroup<4> >
   */
   MPI::Datatype MpiTraits< Simp::SpeciesGroup<4> >::type = MPI::BYTE;
   bool MpiTraits< Simp::SpeciesGroup<4> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<3> >
   */
   MPI::Datatype MpiTraits< Simp::SpeciesGroup<3> >::type = MPI::BYTE;
   bool MpiTraits< Simp::SpeciesGroup<3> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<2> >
   */
   MPI::Datatype MpiTraits< Simp::SpeciesGroup<2> >::type = MPI::BYTE;
   bool MpiTraits< Simp::SpeciesGroup<2> >::hasType = false;

}

namespace Simp
{

   void commitMpiSpeciesGroupTypes()
   {
      if (!Util::MpiTraits< SpeciesGroup<2> >::hasType) {
         SpeciesGroup<2>::commitMpiType();
      }
      if (!Util::MpiTraits< SpeciesGroup<3> >::hasType) {
         SpeciesGroup<3>::commitMpiType();
      }
      if (!Util::MpiTraits< SpeciesGroup<4> >::hasType) {
         SpeciesGroup<4>::commitMpiType();
      }
   }

}
#endif // ifdef  UTIL_MPI
