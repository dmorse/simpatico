#ifndef MCMD_MASK_POLICY_H
#define MCMD_MASK_POLICY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   /**
   * Enumeration of policies for suppressing ("masking") some pair interactions
   *
   * Possible values:
   *
   *  - MaskNone:   all atoms can undergo pair interactions.
   *  - MaskBonded: mask nonbonded pair interactions between bonded atoms.
   *
   * \ingroup McMd_Chemistry_Module
   */
   enum MaskPolicy {MaskNone=0, MaskBonded=1};

   /**
   * istream extractor for a MaskPolicy.
   *
   * \param  in       input stream
   * \param  policy   MaskPolicy to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, MaskPolicy& policy);

   /**
   * ostream inserter for an MaskPolicy.
   *
   * \param  out      output stream
   * \param  policy   MaskPolicy to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, MaskPolicy policy);

   /**
   * Serialize a MaskPolicy
   *
   * \param ar      archive object
   * \param policy  object to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, MaskPolicy& policy, const unsigned int version)
   { serializeEnum(ar, policy, version); }

}

#include <util/global.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Explicit specialization MpiTraits<MaskPolicy>.
   */
   template <>
   class MpiTraits<McMd::MaskPolicy> {  
   public:  
      static MPI::Datatype type;     ///< MPI Datatype
      static bool hasType;           ///< Is the MPI type initialized?
   };

}
#endif

#endif
