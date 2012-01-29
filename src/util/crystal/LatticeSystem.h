#ifndef LATTICE_SYSTEM_H
#define LATTICE_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif

#include <iostream>

namespace Util
{

   /**
   * Enumeration of the 7 possible Bravais lattice systems.
   *
   * \ingroup Crystal_Module
   */
   enum LatticeSystem {Cubic, Tetragonal, Orthorhombic, 
                       Monoclinic, Triclinic, Rhombohedral, Hexagonal};


   /**
   * istream extractor for a LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, LatticeSystem& lattice);

   /**
   * ostream inserter for an LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, LatticeSystem lattice);

   #ifdef UTIL_MPI

   /**
   * Explicit specialization MpiTraits<LatticeSystem>.
   */
   template <>
   class MpiTraits<Util::LatticeSystem>
   {  
   public:  
      static MPI::Datatype type;      ///< MPI Datatype
      static bool hasType;            ///< Is the MPI type initialized?
   };

   #endif
}
#endif
