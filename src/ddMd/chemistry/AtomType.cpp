#ifndef DDMD_ATOM_TYPE_CPP
#define DDMD_ATOM_TYPE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomType.h"
#include <util/param/Parameter.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor
   */
   AtomType::AtomType()
    : mass_(1.0),
      name_(),
      id_(-1) 
   {}

   /*
   * Set the index.
   */
   void AtomType::setId(int id)
   {  id_ = id; }

   /* 
   * Input a AtomType from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, AtomType &atomType)
   {
      in >> atomType.name_;
      in >> atomType.mass_;
      return in;
   }
   
   /* 
   * Output a AtomType to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const AtomType &atomType) 
   {
      out.width(Parameter::Width);
      out << atomType.name_;
      out.setf(std::ios::scientific);
      out.width(Parameter::Width);
      out.precision(Parameter::Precision);
      out << atomType.mass_;
      return out;
   }

   #ifdef UTIL_MPI
   /*
   * Call to guarantee initialization of static data.
   */
   void AtomType::initStatic()
   {
      MpiTraits<DdMd::AtomType>::type = MPI::BYTE;
      MpiTraits<DdMd::AtomType>::hasType = false;
   }
   #endif

} 

#ifdef UTIL_MPI
namespace Util
{

   template <>
   void send<DdMd::AtomType>(MPI::Comm& comm, DdMd::AtomType& data, int dest, int tag)
   {
      std::string name = data.name();
      double      mass = data.mass();
      send<std::string>(comm, name, dest, tag);
      send<double>(comm, mass, dest, tag);
   }

   template <>
   void recv<DdMd::AtomType>(MPI::Comm& comm, DdMd::AtomType& data, int source, int tag)
   {
      std::string name;
      double      mass; 
      recv<std::string>(comm, name, source, tag);
      recv<double>(comm, mass, source, tag);
      data.setName(name);
      data.setMass(mass);
   }

   template <>
   void bcast<DdMd::AtomType>(MPI::Intracomm& comm, DdMd::AtomType& data, int root)
   {
      std::string name;
      double      mass; 
      int         rank = comm.Get_rank();
      if (rank == root) {
         name = data.name();
         mass = data.mass();
      }
      bcast<std::string>(comm, name, root);
      bcast<double>(comm, mass, root);
      if (rank != root) {
         data.setName(name);
         data.setMass(mass);
      }
   }

   /**
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<DdMd::AtomType>::type = MPI::BYTE;
   bool MpiTraits<DdMd::AtomType>::hasType = false;

}
#endif

#endif
