#ifndef MCMD_MPI_H
#define MCMD_MPI_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#ifdef UTIL_MPI

#include <util/mpi/MpiTraits.h>
#include <util/containers/Pair.h>

namespace Util
{

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
