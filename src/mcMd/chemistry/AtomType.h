#ifndef MCMD_ATOM_TYPE_H
#define MCMD_ATOM_TYPE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <string>
#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * Descriptor for a type of Atom.
   *
   * An AtomTypeType has a mass, a name string, and an integer id.
   * If coulomb interactions are enabled (ifdef SIMP_COULOMB), it
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
      * Set the name string.
      *
      * \param name name string
      */
      void setName(std::string name);

      /**
      * Set the mass.
      *
      * \param mass atom mass
      */
      void setMass(double mass);

      #ifdef SIMP_COULOMB
      /**
      * Set the boolean "hasCharge" property.
      *
      * An AtomType has an associated electrical charge value if and
      * only if hasCharge is true.  A charge value appears in the
      * text file format used by the inserter and extractor iostream
      * operators if and only if the hasCharge property is set true. 
      *
      * The hasCharge property should be set true for all atom types 
      * (even those with no charge) if the system has any charged
      * atom types, and thus has Coulomb interactions, and should be
      * set false for all atom types for a neutral system with no 
      * Coulomb interactions. 
      *
      * \param hasCharge true if this system has Coulomb interactions.
      */
      void setHasCharge(bool hasCharge);

      /**
      * Set the charge value.
      *
      * Precondition: The hasCharge property must have been set true.
      *
      * \param charge atom electrical charge
      */
      void setCharge(double charge);
      #endif

      //@}
      /// \name Accessors
      //@{

      /**
      * Get the mass.
      */
      double mass() const;

      #ifdef SIMP_COULOMB
      /**
      * Does this type have a charge value?
      */
      bool hasCharge() const;

      /**
      * Get the electrical charge value.
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

      #ifdef SIMP_COULOMB
      /// Electrical charge of atom type.
      double  charge_;
      #endif

      /// Name of type.
      std::string  name_;

      /// Integer index.
      int  id_;

      /// Does this type have a charge value?
      bool hasCharge_;

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
   * Get the type id.
   */
   inline int  AtomType::id() const
   {  return id_; }

   /*
   * Get the name string.
   */
   inline const std::string& AtomType::name() const
   {  return name_; }

   /*
   * Get the mass.
   */
   inline double AtomType::mass() const
   {  return mass_; }

   #ifdef SIMP_COULOMB
   /*
   * Get the hasCharge bool flag.
   */
   inline bool AtomType::hasCharge() const
   {  return hasCharge_; }

   /*
   * Get the electrical charge.
   */
   inline double AtomType::charge() const
   {
      UTIL_ASSERT(hasCharge_);  
      return charge_; 
   }
   #endif

   // Declarations of friend extracter and inserter functions

   /**
   * istream extractor (>>) for an AtomType.
   *
   * Format:
   *
   *    name  [string] mass [double] charge [double]
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
   *    name  mass [charge]
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
      #ifdef SIMP_COULOMB
      ar & atomType.hasCharge_;
      if (atomType.hasCharge_) {
         ar & atomType.charge_;
      }
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
