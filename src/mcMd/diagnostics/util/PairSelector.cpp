#ifndef PAIR_SELECTOR_CPP
#define MCMD_PAIR_SELECTOR_CPP

#include "PairSelector.h"
#include <mcMd/chemistry/Atom.h>
#include <util/param/Parameter.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   PairSelector::PairSelector()
    : pairType_(ALL),
      atom1TypeId_(-1),
      atom2TypeId_(-1),
      avoidDoubleCounting_(false)
   {}

   /*
   * Set policy to avoid double counting (true) or not (true).
   */
   void PairSelector::setAvoidDoubleCounting(bool avoidDoubleCounting)
   {  avoidDoubleCounting_ = avoidDoubleCounting; }

   /**
   * Return true to accept this atom pair, false to reject.
   */
   bool PairSelector::match(const Atom& atom1, const Atom& atom2) const
   {

      // Check molecular identities
      if (pairType_ == INTRA && &atom1.molecule() != &atom2.molecule() ) {
         return false;
      } else
      if (pairType_ == INTER && &atom1.molecule() == &atom2.molecule() ) {
         return false;
      }

      // Exclude case where atom1 and atom2 are identical
      if (&atom1 == &atom2) {
         return false;
      }

      // Check atom types
      if (atom1TypeId_ >= 0 && atom1.typeId() != atom1TypeId_) {
         return false;
      }
      if (atom2TypeId_ >= 0 && atom2.typeId() != atom2TypeId_) {
         return false;
      }

      // Check for double counting
      if (avoidDoubleCounting_) {
         if (atom1TypeId_ < 0 && atom2TypeId_ < 0) {
            if (atom1.id() > atom2.id()) {
               return false;
            }		  
         } else 
         if (atom1.typeId() == atom2.typeId()) {
            if (atom1.id() > atom2.id()) {
               return false;
            }		  
         }
      }
     
      // If the pair has passed all previous tests, accept it. 
      return true;

   }

   /* 
   * Input a PairSelector from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, PairSelector &selector)
   {
      in >> selector.pairType_;
      in >> selector.atom1TypeId_;
      in >> selector.atom2TypeId_;
      return in;
   }
   
   /* 
   * Output a PairSelector to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const PairSelector &selector) 
   {
      out.width(Parameter::Width);
      out << selector.pairType_;
      out.width(10);
      out << selector.atom1TypeId_;
      out.width(10);
      out << selector.atom2TypeId_;
      return out;
   }

   /* 
   * Input a PairSelector::PairType from an istream.
   */
   std::istream& operator>>(std::istream& in, PairSelector::PairType &type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "intra" || buffer == "Intra" || buffer == "INTRA") {
         type = PairSelector::INTRA;
      } else
      if (buffer == "inter" || buffer == "Inter" || buffer == "INTER") {
         type = PairSelector::INTER;
      } else
      if (buffer == "all" || buffer == "All" || buffer == "ALL") {
         type = PairSelector::ALL;
      } else {
         UTIL_THROW("Unknown PairSelector::Type");
      }	      
      return in;
   }
   
   /* 
   * Output a PairSelector to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const PairSelector::PairType &type) 
   {
      if (type == PairSelector::INTRA) {
         out << "Intra"; 
      } else
      if (type == PairSelector::INTER) {
         out << "Inter"; 
      } else
      if (type == PairSelector::ALL) {
         out << "All"; 
      } else {
         UTIL_THROW("Unknown PairSelector::PairType");
      }	      
      return out;
   }

} 

#ifdef UTIL_MPI

namespace Util
{

   // Initialize MpiTraits<McMd::PairSelector>
   MPI::Datatype MpiTraits< McMd::PairSelector>::type = MPI::BYTE;
   bool MpiTraits< McMd::PairSelector>::hasType = false;

   // Initialize MpiTraits<McMd::PairSelector::PairType>
   MPI::Datatype MpiTraits< McMd::PairSelector::PairType>::type = MPI::INT;
   bool MpiTraits< McMd::PairSelector::PairType>::hasType = true;

}

#include <util/mpi/MpiStructBuilder.h>   

namespace McMd
{

   /*
   * Commit MPI Datatype.
   */
   void PairSelector::commitMpiType() 
   {
      MpiStructBuilder builder;
      PairSelector     object;
      builder.setBase(&object);
      builder.addMember(&object.pairType_, MpiTraits<PairSelector::PairType>::type);
      builder.addMember(&object.atom1TypeId_, MPI::INT);
      builder.addMember(&object.atom2TypeId_, MPI::INT);
      builder.commit(Util::MpiTraits<PairSelector>::type);
      Util::MpiTraits<PairSelector>::hasType = true;
   }

}
#endif

#endif
