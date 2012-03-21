#ifndef DDMD_BOUNDARY_ENSEMBLE_CPP
#define DDMD_BOUNDARY_ENSEMBLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BoundaryEnsemble.h"

#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>

namespace Util{

   /**
   * Initialize BoundaryEnsemble MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::BoundaryEnsemble>::type = MPI::BYTE;
   bool MpiTraits<DdMd::BoundaryEnsemble>::hasType = false;

   /**
   * Initialize BoundaryEnsemble::Type MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::BoundaryEnsemble::Type>::type = MPI::INT;
   bool MpiTraits<DdMd::BoundaryEnsemble::Type>::hasType = true;

}

#endif // ifdef UTIL_MPI

namespace DdMd
{

   using namespace Util;

   /**
   * Constructor.
   */
   BoundaryEnsemble::BoundaryEnsemble(Type type)
    : pressure_(1.0),
      type_(type)
   {}

   /**
   * Set the pressure.
   */
   void  BoundaryEnsemble::setPressure(double pressure)
   {
      if (!isIsobaric()) {
	 UTIL_THROW("Must be an isobaric ensemble");
      }
      pressure_ = pressure;
   }

   /**
   * Read the type and (if necessary) pressure from file.
   */
   void BoundaryEnsemble::readParam(std::istream& in)
   {
      readBegin(in, "BoundaryEnsemble");
      read<Type>(in, "type", type_);
      if (isIsobaric()) {
         read<double>(in, "pressure", pressure_);
      }
      readEnd(in);
   }

   /*
   * Extract an BoundaryEnsemble::Type from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, BoundaryEnsemble::Type &type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "RIGID" || buffer == "rigid") {
         type = BoundaryEnsemble::RIGID;
      } else
      if (buffer == "ISOBARIC" || buffer == "isobaric") {
         type = BoundaryEnsemble::ISOBARIC;
      } else {
         UTIL_THROW("Invalid BoundaryEnsemble::Type value input");
      }
      return in;
   }

   /*
   * Insert a BoundaryEnsemble::Type to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, const BoundaryEnsemble::Type &type)
   {
      if (type == BoundaryEnsemble::RIGID) {
         out << "rigid";
      } else
      if (type == BoundaryEnsemble::ISOBARIC) {
         out << "isobaric";
      } else
      if (type == BoundaryEnsemble::UNKNOWN) {
         out << "unknown";
      }
      return out;
   }

   #ifdef UTIL_MPI

   /**
   * Commit MPI Datatype.
   */
   void BoundaryEnsemble::commitMpiType()
   {
      MpiStructBuilder builder;
      BoundaryEnsemble   object;

      builder.setBase(&object);
      builder.addMember(&object.pressure_, MPI::DOUBLE);
      builder.addMember(&object.type_, MPI::INT);
      builder.commit(Util::MpiTraits<BoundaryEnsemble>::type);
      Util::MpiTraits<DdMd::BoundaryEnsemble>::hasType = true;
   }

   #endif

} //namespace DdMd
#endif
