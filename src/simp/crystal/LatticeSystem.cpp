/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LatticeSystem.h"      // class header
#include <util/misc/ioUtil.h> 
#include <util/global.h>      

namespace Simp
{

   using namespace Util; 

   #ifdef UTIL_MPI
   /**
   * Initialize LatticeSystem:: MPI Datatype.
   */
   MPI::Datatype MpiTraits<Simp::LatticeSystem>::type = MPI::INT;
   bool MpiTraits<Simp::LatticeSystem>::hasType = true;
   #endif

   /* 
   * Extract a LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, LatticeSystem& lattice)
   {
      // Precondition (check state of input stream)
      try {
         checkRequiredIstream(in);
      } catch (Exception e) {
         Log::file() << "Error in istream reading LatticeSystem" << std::endl;
         throw e;
      }

      // Read string representation
      std::string buffer;
      skipws(in);
      in >> buffer;
      if (buffer.size() == 0) {
         UTIL_THROW("Empty LatticeSystem string after read");
      }

      // Match string representation and set value
      if (buffer == "Cubic" || buffer == "cubic") {
         lattice = Cubic;
      } else 
      if (buffer == "Tetragonal" || buffer == "tetragonal") {
         lattice = Tetragonal;
      } else 
      if (buffer == "Orthorhombic" || buffer == "orthorhombic") {
         lattice = Orthorhombic;
      } else 
      if (buffer == "Monoclinic" || buffer == "monoclinic") {
         lattice = Monoclinic; 
      } else 
      if (buffer == "Triclinic" || buffer == "triclinic") {
         lattice = Triclinic; 
      } else 
      if (buffer == "Rhombohedral" || buffer == "rhombohedral") {
         lattice = Rhombohedral;
      } else 
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = Hexagonal; 
      } else {
         #ifndef UTIL_MPI
         Log::file() << "Unknown LatticeSystem. Value =|" 
                     << buffer << "|" << std::endl;
         #endif
         UTIL_THROW("Invalid LatticeSystem value input");
      }
      return in;
   }
   
   /* 
   * Insert a LatticeSystem to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, LatticeSystem lattice) 
   {
      if (lattice == Cubic) {
         out << "cubic";
      } else 
      if (lattice == Tetragonal) {
         out << "tetragonal";
      } else
      if (lattice == Orthorhombic) {
         out << "orthorhombic";
      } else
      if (lattice == Monoclinic) {
         out << "monoclinic";
      } else
      if (lattice == Triclinic) {
         out << "triclinic";
      } else
      if (lattice == Rhombohedral) {
         out << "rhombohedral";
      } else
      if (lattice == Hexagonal) {
         out << "hexagonal";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out; 
   }

}
