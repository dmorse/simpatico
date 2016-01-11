/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairSelector.h"
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomContext.h>
#include <util/param/Parameter.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif


namespace DdMd
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
      if (pairType_ != ALL) {
         if (not (atom1.hasAtomContext() && atom2.hasAtomContext())) {
            UTIL_THROW("Pairtype != ALL and hasAtomContext not enabled");
         }

         // Check molecular identities
         AtomContext ac1 = atom1.context();
         AtomContext ac2 = atom2.context();

         bool sameMolecule = (atom1.context().speciesId  == atom2.context().speciesId);
         sameMolecule = sameMolecule &&
                  (atom1.context().moleculeId == atom2.context().moleculeId);

         if (pairType_ == INTRA && not(sameMolecule)) {
            return false;
         }
         if (pairType_ == INTER && sameMolecule) {
            return false;
         }
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

   #ifdef UTIL_MPI
   /*
   * Call to guarantee initialization of static data.
   */
//   void PairSelector::initStatic()
//   {
//      MpiTraits<DdMd::PairSelector>::type = MPI::BYTE;
//      MpiTraits<DdMd::PairSelector>::hasType = false;
//   }
   #endif

}

   #ifdef UTIL_MPI
namespace Util
{

   template <>
   void send<DdMd::PairSelector>(MPI::Comm& comm, DdMd::PairSelector& data, int dest, int tag)
   {
      DdMd::PairSelector::PairType         pairType = data.pairType();
      int      atom1TypeId = data.atom1TypeId();
      int      atom2TypeId = data.atom2TypeId();
      bool           avoid = data.avoidDoubleCounting();

      send<DdMd::PairSelector::PairType>(comm, pairType, dest, tag);
      send<int>(comm, atom1TypeId, dest, tag);
      send<int>(comm, atom2TypeId, dest, tag);
      send<bool>(comm, avoid, dest, tag);
   }

   template <>
   void recv<DdMd::PairSelector>(MPI::Comm& comm, DdMd::PairSelector& data, int source, int tag)
   {
      DdMd::PairSelector::PairType         pairType;
      int         atom1TypeId;
      int         atom2TypeId;
      bool        avoid;
//      std::string name;
//      double      mass; 
      recv<DdMd::PairSelector::PairType>(comm, pairType, source, tag);
      recv<int>(comm, atom1TypeId, source, tag);
      recv<int>(comm, atom2TypeId, source, tag);
      recv<bool>(comm, avoid, source, tag);

      int         rank = comm.Get_rank();
         std::string str;
         switch (pairType) {
            case 1:
               str = "INTRA"; break;
            case 2:
               str = "INTER"; break;
            case 3:
               str = "ALL"; break;
            default:
               str = "shit";
         }
         printf("%i : RECV    type1[%i] type2[%i] pt[%s]\n", rank,
               atom1TypeId, atom2TypeId, str.c_str());


      data.setPairType(pairType);
      data.setAtom1TypeId(atom1TypeId);
      data.setAtom2TypeId(atom2TypeId);
      data.setAvoidDoubleCounting(avoid);
   }
//
   template <>
   void bcast<DdMd::PairSelector>(MPI::Intracomm& comm, DdMd::PairSelector& data, int root)
   {
      DdMd::PairSelector::PairType         pairType;
      int         atom1TypeId;
      int         atom2TypeId;
      bool        avoid;

//      std::string name;
//      double      mass; 
      int         rank = comm.Get_rank();
      if (rank == root) {
         pairType    = data.pairType();
         atom1TypeId = data.atom1TypeId();
         atom2TypeId = data.atom2TypeId();
         avoid       = data.avoidDoubleCounting();
      }
      bcast<DdMd::PairSelector::PairType>(comm, pairType, root);
      bcast<int>(comm, atom1TypeId, root);
      bcast<int>(comm, atom2TypeId, root);
      bcast<bool>(comm, avoid, root);

//      bcast<std::string>(comm, name, root);
//      bcast<double>(comm, mass, root);
      if (rank != root) {
         data.setPairType(pairType);
         data.setAtom1TypeId(atom1TypeId);
         data.setAtom2TypeId(atom2TypeId);
         data.setAvoidDoubleCounting(avoid);
//         data.setName(name);
//         data.setMass(mass);
      }
   }

   /**
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::PairSelector>::type = MPI::BYTE;
   bool MpiTraits<DdMd::PairSelector>::hasType = false;


   MPI::Datatype MpiTraits<DdMd::PairSelector::PairType>::type = MPI::INT;
   bool MpiTraits<DdMd::PairSelector::PairType>::hasType = true;

}


   /*
   * Commit MPI Datatype.
   */
//   void PairSelector::commitMpiType() 
//   {
//      if (!MpiTraits<PairSelector>::hasType) {
//         MpiStructBuilder builder;
//         PairSelector selector;
//         builder.setBase(&selector);
//         builder.addMember(&atom1TypeId_, MPI::INT);
//         builder.addMember(&atom2TypeId_, MPI::INT);
//         builder.addMember(&pairType_, MPI::INT);
//         builder.addMember(&avoidDoubleCounting_, MPI::BOOL);
//         builder.commit(MpiTraits<PairSelector>::type);
//         MpiTraits<PairSelector>::hasType = true;
//      }
//   }
   #endif


