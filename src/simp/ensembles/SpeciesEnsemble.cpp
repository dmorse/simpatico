/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesEnsemble.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpeciesEnsemble::SpeciesEnsemble(Type type)
    : mu_(1.0),
      type_(type)
   {}

   /*
   * Destructor.
   */
   SpeciesEnsemble::~SpeciesEnsemble()
   {}

   /*
   * Set the chemical potential mu.
   */
   void  SpeciesEnsemble::setMu(double mu)
   { 
      if (!isGrand()) {
         UTIL_THROW("Must be a grand-canonical ensemble");
      }
      mu_ = mu; 
   }

   /*
   * Read the type and (if necessary) chemical potential mu from file.
   */
   void SpeciesEnsemble::readParam(std::istream& in)
   {
      readBegin(in, "SpeciesEnsemble");
      read<Type>(in, "type", type_);
      if (isGrand()) {
         read<double>(in, "mu", mu_); 
      }
      readEnd(in);
   } 

   /* 
   * Extract an SpeciesEnsemble::Type from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, SpeciesEnsemble::Type &type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "CLOSED" || buffer == "closed") {
         type = SpeciesEnsemble::CLOSED;
      } else 
      if (buffer == "GRAND" || buffer == "grand") {
         type = SpeciesEnsemble::GRAND;
      } else {
         UTIL_THROW("Invalid SpeciesEnsemble::Type value input");
      }
      return in;
   }
   
   /* 
   * Insert a SpeciesEnsemble::Type to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, const SpeciesEnsemble::Type &type) 
   {
      if (type == SpeciesEnsemble::CLOSED) {
         out << "closed";
      } else 
      if (type == SpeciesEnsemble::GRAND) {
         out << "grand";
      }
      return out; 
   }

   #ifdef UTIL_MPI
   /*
   * Commit MPI Datatype.
   */
   void SpeciesEnsemble::commitMpiType() 
   {
      MpiStructBuilder builder;
      SpeciesEnsemble  object;

      builder.setBase(&object);
      builder.addMember(&object.mu_, MPI_DOUBLE);
      builder.addMember(&object.type_, MPI_INT);
      builder.commit(MpiTraits<SpeciesEnsemble>::type);
      MpiTraits<SpeciesEnsemble>::hasType = true;
   }
   #endif

}
 
#ifdef UTIL_MPI
namespace Util
{

   /*
   * Initialize MPI Datatype.
   */
   MPI_Datatype MpiTraits<Simp::SpeciesEnsemble>::type = MPI_BYTE;
   bool MpiTraits<Simp::SpeciesEnsemble>::hasType = false;

   /*
   * Initialize MPI Datatype.
   */
   MPI_Datatype MpiTraits<Simp::SpeciesEnsemble::Type>::type = MPI_INT;
   bool MpiTraits<Simp::SpeciesEnsemble::Type>::hasType = true;

} 
#endif

