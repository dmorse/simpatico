#ifndef MCMD_MASK_POLICY_CPP
#define MCMD_MASK_POLICY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskPolicy.h"    // class header

#ifdef UTIL_MPI
namespace Util
{

   /**
   * Initialize MPI Datatype associated with MaskPolicy.
   */
   MPI::Datatype MpiTraits<McMd::MaskPolicy>::type    = MPI::INT;
   bool          MpiTraits<McMd::MaskPolicy>::hasType = true;

}
#endif

namespace McMd
{

   using namespace Util;

   /* 
   * Extract a MaskPolicy from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, MaskPolicy& policy)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "MaskNone" || buffer == "maskNone") {
         policy = MaskNone;
      } else 
      if (buffer == "MaskBonded" || buffer == "maskBonded") {
         policy = MaskBonded;
      } else {
         UTIL_THROW("Invalid MaskPolicy string");
      } 
      return in;
   }
   
   /* 
   * Insert a MaskPolicy to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, MaskPolicy policy) 
   {
      if (policy == MaskNone) {
         out << "MaskNone";
      } else 
      if (policy == MaskBonded) {
         out << "MaskBonded";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out; 
   }

}
#endif
