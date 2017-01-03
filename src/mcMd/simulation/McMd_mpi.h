#ifndef MCMD_MPI_H
#define MCMD_MPI_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#ifdef UTIL_MPI

#include <util/mpi/MpiTraits.h>
#include <mcMd/species/SpeciesGroup.h>
#include <util/containers/Pair.h>

namespace Util
{
   /**
   * Explicit specialization MpiTraits< SpeciesGroup<4> >.
   */
   template <>
   class MpiTraits< McMd::SpeciesGroup<4> >
   {  
   public:  
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };


   /**
   * Explicit specialization MpiTraits< SpeciesGroup<3> >.
   */
   template <>
   class MpiTraits< McMd::SpeciesGroup<3> >
   {  
   public:  
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };


   /**
   * Explicit specialization MpiTraits< SpeciesGroup<2> >.
   */
   template <>
   class MpiTraits< McMd::SpeciesGroup<2> >
   {  
   public:  
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * Explicit specialization MpiTraits< Pair<int> >.
   */
   template <>
   class MpiTraits< Pair<int> >
   {  
   public:  
      static MPI::Datatype type;  ///< MPI Datatype
      static bool hasType;        ///< Is the MPI type initialized?
   };

} 

namespace McMd
{

   /**
   * Commit all MPI data types needed for Mc and Md simulations.
   */  
   void commitMpiTypes(); 

}
 

#endif  // ifdef UTIL_MPI
#endif  // ifndef MCMD_MPI_H
