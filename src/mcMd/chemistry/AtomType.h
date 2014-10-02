#ifndef MCMD_ATOM_TYPE_H
#define MCMD_ATOM_TYPE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <string>
#include <iostream>

namespace McMd
{

   /**
   * Descriptor for a type of Atom.
   *
   * An AtomTypeType has a mass, a name string, and an integer id.
   * If coulomb interactions are enabled (ifdef INTER_COULOMB), it
   * also has a electrical charge.
   *
   * \ingroup McMd_Chemistry_Module
   */
   class AtomType
   {

   public:

      /**
      * Constructor.
      */
      AtomType();

      /**
      * Destructor.
      */
      virtual ~AtomType();

      /// \name Mutators
      //@{

      /**
      * Set the type index.
      *
      * \param Id integer index.
      */
      void setId(int Id);

      /**
      * Set the mass.
      *
      * \param mass atom mass
      */
      void setMass(double mass);

      #ifdef INTER_COULOMB
      /**
      * Set the charge.
      *
      * \param charge atom electrical charge
      */
      void setCharge(double charge);
      #endif

      /**
      * Set the name string.
      *
      * \param name name string
      */
      void setName(std::string name);

      //@}
      /// \name Accessors
      //@{

      /**
      * Get the mass.
      */
      double mass() const;

      #ifdef INTER_COULOMB
      /**
      * Get the electrical charge.
      */
      double charge() const;
      #endif

      /**
      * Get the name string.
      */
      const std::string& name() const;

      /// Get the index.
      int id() const;

      //@}

   private:

      /// Mass of atom type.
      double  mass_;

      #ifdef INTER_COULOMB
      /// Electrical charge of atom type.
      double  charge_;
      #endif

      /// Name of type.
      std::string  name_;

      /// Integer index.
      int  id_;

   //friends:

      friend 
      std::istream& operator>>(std::istream& in, AtomType &atomType);

      friend 
      std::ostream& operator<<(std::ostream& out, const AtomType &atomType);

      template <class Archive> 
      friend void 
      serialize(Archive& ar, AtomType& atomType, const unsigned int version);

   };

   // Inline member functions.

   /*
   * Set the mass.
   */
   inline void AtomType::setMass(double mass)
   {  mass_ = mass; }

   #ifdef INTER_COULOMB
   /*
   * Set the electrical charge.
   */
   inline void AtomType::setCharge(double charge)
   {  charge_ = charge; }
   #endif 

   /*
   * Set the name string.
   */
   inline void AtomType::setName(std::string name)
   {  name_ = name; }

   /*
   * Get the mass.
   */
   inline double AtomType::mass() const
   {  return mass_; }

   #ifdef INTER_COULOMB
   /*
   * Get the electrical charge.
   */
   inline double AtomType::charge() const
   {  return charge_; }
   #endif

   /*
   * Get the name string.
   */
   inline const std::string& AtomType::name() const
   {  return name_; }

   /*
   * Get the type id.
   */
   inline int  AtomType::id() const
   {  return id_; }

   // Friend operator declarations

   /**
   * istream extractor (>>) for an AtomType.
   *
   * Format:
   *
   *    name  [string] mass [double]
   *
   * \param in  input stream
   * \param atomType  AtomType to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, AtomType &atomType);

   /**
   * ostream inserter (<<) for an AtomType.
   *
   * Format, one one line with no line break:
   *
   *    name  mass
   *
   * \param  out  output stream
   * \param  atomType  AtomType to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const AtomType &atomType);

   /**
   * Serialize an AtomType.
   *
   * \param ar  archive object
   * \param atomType  object to be serialized
   * \param version  archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, AtomType& atomType, const unsigned int version)
   {
      ar & atomType.name_;
      ar & atomType.mass_;
      #ifdef INTER_COULOMB
      ar & atomType.charge_;
      #endif
   }

}

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Send an AtomType (wrapper for MPI Send).
   *
   * \param comm MPI communicator
   * \param data AtomType data
   * \param dest MPI rank of destination (receiving) processor
   * \param tag  integer identifier for message
   */
   template <>
   void send<McMd::AtomType>(MPI::Comm& comm, McMd::AtomType& data, int dest, int tag);

   /**
   * Receive an AtomType (wrapper for MPI Recv).
   *
   * \param comm MPI communicator
   * \param data AtomType data
   * \param source MPI rank of source (sending) processor
   * \param tag  integer identifier for message
   */
   template <>
   void recv<McMd::AtomType>(MPI::Comm& comm, McMd::AtomType& data, int source, int tag);

   /**
   * Broadcast an AtomType (wrapper for MPI Bcast).
   *
   * \param comm MPI communicator
   * \param data AtomType data
   * \param root MPI rank of root processor from which data is broadcast
   */
   template <>
   void bcast<McMd::AtomType>(MPI::Intracomm& comm, McMd::AtomType& data, int root);

   /**
   * Explicit specialization MpiTraits<AtomType>.
   */
   template <>
   class MpiTraits<McMd::AtomType>
   {
   public:
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };

}
#endif

#endif
