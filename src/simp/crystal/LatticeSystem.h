#ifndef SIMP_LATTICE_SYSTEM_H
#define SIMP_LATTICE_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif

#include <iostream>

namespace Simp
{

   /**
   * Enumeration of the 7 possible Bravais lattice systems.
   *
   * Allowed values are: Cubic, Tetragonal, Orthorhombic, Monoclinic
   * Triclinic, Rhombohedral, and Hexagonal.
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

   /**
   * Serialize a LatticeSystem value.
   *
   * \param ar      archive object
   * \param lattice value to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, LatticeSystem& lattice, const unsigned int version)
   {  serializeEnum(ar, lattice, version); }

}

#ifdef UTIL_MPI
namespace Util 
{

   /**
   * Explicit specialization MpiTraits<LatticeSystem>.
   */
   template <>
   class MpiTraits<Simp::LatticeSystem>
   {
   public:
      static MPI::Datatype type;  ///< MPI Datatype
      static bool hasType;        ///< Is the MPI type initialized?
   };

}
#endif

#endif
