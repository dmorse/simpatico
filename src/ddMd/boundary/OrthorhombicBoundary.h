#ifndef ORTHORHOMBIC_BOUNDARY_H
#define ORTHORHOMBIC_BOUNDARY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoBoundaryBase.h"         // base class
#include <iostream>

class OrthorhombicBoundaryTest;

namespace DdMd
{

   using namespace Util;

   /**
   * An orthorhombic periodic unit cell.
   *
   * \ingroup Boundary_Module
   */
   class OrthorhombicBoundary : public OrthoBoundaryBase
   {

   public:

      /**
      * Constructor.
      */
      OrthorhombicBoundary();

      /**
      * Set unit cell dimensions.
      *
      * Also sets all related lengths and volume.
      *
      * \param lengths  Vector of unit cell lengths. 
      */
      void setLengths(const Vector &lengths);

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

   // friends:

      /// Unit test
      friend class ::OrthorhombicBoundaryTest;

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, OrthorhombicBoundary& boundary);

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, const OrthorhombicBoundary& boundary);

   };

   /**
   * istream extractor for a OrthorhombicBoundary.
   *
   * \param  in        input stream
   * \param  boundary  OrthorhombicBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, OrthorhombicBoundary& boundary);

   /**
   * ostream inserter for an OrthorhombicBoundary.
   *
   * \param  out      output stream
   * \param  boundary OrthorhombicBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const OrthorhombicBoundary& boundary);

}
 
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>
#include <mpi.h>

namespace Util
{

   template <>
   void send<DdMd::OrthorhombicBoundary>(MPI::Comm& comm, 
             DdMd::OrthorhombicBoundary& data, int dest, int tag);

   template <>
   void recv<DdMd::OrthorhombicBoundary>(MPI::Comm& comm, 
             DdMd::OrthorhombicBoundary& data, int source, int tag);

   template <>
   void bcast<DdMd::OrthorhombicBoundary>(MPI::Intracomm& comm, 
              DdMd::OrthorhombicBoundary& data, int root);

   /**
   * Explicit specialization MpiTraits<OrthorhombicBoundary>.
   */
   template <>
   class MpiTraits<DdMd::OrthorhombicBoundary>
   {
   public:
      static MPI::Datatype type;
      static bool hasType;
   };

}
#endif

#endif
