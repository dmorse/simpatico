#ifndef UTIL_BOUNDARY_SERIALIZE_H
#define UTIL_BOUNDARY_SERIALIZE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoRegion.h"  
#include "OrthorhombicBoundary.h"  
#include "MonoclinicBoundary.h"  

namespace Util
{

   /*
   * Serialize an OrthoRegion to/from an archive.
   */
   template <class Archive>
   void OrthoRegion::serialize(Archive& ar, const unsigned int version)
   {
      ar & minima_;
      ar & maxima_;
      ar & lengths_;
      ar & halfLengths_;
      ar & volume_;
   }

   /*
   * Serialize an OrthorhombicBoundary to/from an archive.
   */
   template <class Archive>
   void 
   OrthorhombicBoundary::serialize(Archive& ar, const unsigned int version)
   {
      OrthoRegion::serialize(ar, version);
      serializeEnum(ar, lattice_, version);
      if (Archive::is_loading()) {
         reset();
         isValid();
      }
   }

   /*
   * Serialize an OrthorhombicBoundary to/from an archive.
   */
   template <class Archive>
   void 
   MonoclinicBoundary::serialize(Archive& ar, const unsigned int version)
   {
      ar & l_;
      ar & d_;
      if (Archive::is_loading()) {
         reset();
         isValid();
      }
   }

}
#endif
