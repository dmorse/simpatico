#ifndef SIMP_SPECIES_GROUP_TPP
#define SIMP_SPECIES_GROUP_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesGroup.h"
#include <util/param/Parameter.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>   
#endif

namespace Simp
{

   using namespace Util;

   /* 
   * Input a SpeciesGroup<N> from an istream, without line breaks.
   */
   template <int NAtom>
   std::istream& operator >> (std::istream& in, SpeciesGroup<NAtom>& speciesGroup)
   {
      for (int i = 0; i < NAtom; ++i) {
         in >> speciesGroup.atomIds_[i];
      }
      in >> speciesGroup.typeId_;
      return in;
   }
   
   /* 
   * Output a SpeciesGroup<N> to an ostream, without line breaks.
   */
   template <int NAtom>
   std::ostream& 
   operator << (std::ostream& out, const SpeciesGroup<NAtom>& speciesGroup) 
   {
      for (int i = 0; i < NAtom; ++i) {
         if (i == 0) {
            out.width(Parameter::Width);
         } else {
            out.width(10);
         }
         out << speciesGroup.atomIds_[i];
      }
      out.width(10);
      out << speciesGroup.typeId_;
      return out;
   }

   #ifdef UTIL_MPI
   /**
   * Commit MPI Datatype.
   */
   template <int NAtom>
   void SpeciesGroup<NAtom>::commitMpiType() 
   {
      if (!Util::MpiTraits< SpeciesGroup<NAtom> >::hasType) {
         MpiStructBuilder    builder;
         SpeciesGroup<NAtom> object;
   
         builder.setBase(&object);
         for (int i = 0; i < NAtom; ++i) {
            builder.addMember(&(object.atomIds_[i]), MPI::INT);
         }
         builder.addMember(&object.typeId_, MPI::INT);
         builder.commit(Util::MpiTraits< SpeciesGroup<NAtom> >::type);
         Util::MpiTraits< SpeciesGroup<NAtom> >::hasType = true;
      }
   }
   #endif

}
#endif
