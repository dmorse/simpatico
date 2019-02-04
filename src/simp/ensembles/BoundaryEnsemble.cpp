/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoundaryEnsemble.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   BoundaryEnsemble::BoundaryEnsemble(Type type)
    : pressure_(1.0),
      type_(type)
   {  setClassName("BoundaryEnsemble"); }

   /*
   * Destructor.
   */
   BoundaryEnsemble::~BoundaryEnsemble()
   {}

   /*
   * Set the target pressure.
   */
   void BoundaryEnsemble::setPressure(double pressure)
   {
      if (!isIsobaric()) {
	 UTIL_THROW("Must be an isobaric ensemble");
      }
      pressure_ = pressure;
   }

   /*
   * Read the type and (if necessary) pressure from file.
   */
   void BoundaryEnsemble::readParameters(std::istream& in)
   {
      read<Type>(in, "type", type_);
      if (isIsobaric()) {
         read<double>(in, "pressure", pressure_);
      }
   }

   /*
   * Load internal state from an archive.
   */
   void BoundaryEnsemble::loadParameters(Serializable::IArchive &ar)
   { 
      loadParameter<Type>(ar, "type", type_);
      if (isIsobaric()) {
         loadParameter<double>(ar, "pressure", pressure_);
      }
   }

   /*
   * Save internal state to an archive.
   */
   void BoundaryEnsemble::save(Serializable::OArchive &ar)
   { 
      ar << type_;
      if (isIsobaric()) {
         ar << pressure_;
      }
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
      builder.addMember(&object.pressure_, MPI_DOUBLE);
      builder.addMember(&object.type_, MPI_INT);
      builder.commit(MpiTraits<BoundaryEnsemble>::type);
      MpiTraits<BoundaryEnsemble>::hasType = true;
   }
   #endif

}

#ifdef UTIL_MPI
namespace Util
{

   // Initialize Simp::BoundaryEnsemble MPI Datatype.
   MPI_Datatype MpiTraits<Simp::BoundaryEnsemble>::type = MPI_BYTE;
   bool MpiTraits<Simp::BoundaryEnsemble>::hasType = false;

   // Initialize Simp::BoundaryEnsemble::Type MPI Datatype.
   MPI_Datatype MpiTraits<Simp::BoundaryEnsemble::Type>::type = MPI_INT;
   bool MpiTraits<Simp::BoundaryEnsemble::Type>::hasType = true;

}
#endif
