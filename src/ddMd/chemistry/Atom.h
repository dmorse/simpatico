#ifndef DDMD_ATOM_H
#define DDMD_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#define UTIL_32BIT

#include <util/space/Vector.h>            // members
#include <ddMd/chemistry/Mask.h>          // member
#include <ddMd/communicate/Plan.h>        // member 
#include "AtomArray.h"                    // inline methods

namespace Util{ class Memory; }

namespace DdMd
{

   class Buffer;
   using namespace Util;

   /**
   * A point particle in an MD simulation.
   *
   * An Atom has:
   *
   *   - a position Util::Vector
   *   - a force Util::Vector
   *   - a velocity Util::Vector
   *   - an integer atom type Id
   *   - a boolean isGhost flag
   *   - a global integer id
   *   - a Mask (list of other atoms with masked pair interactions)
   *   - a communication Plan
   *
   * An Atom may only be constructed as an element of an AtomArray. The
   * Atom constructor is private, but is accessible by the AtomArray class
   * via a friend declaration.
   *
   * The interface of the Atom class provides access to each Atom data
   * field through an accessor function, as if all were true C++ class 
   * member. In fact, some of the above are "pseudo-members" that are 
   * stored in separate arrays.  All arrays that store atom data are
   * private members of the associated AtomArray, all of which use the 
   * same indexing scheme to identify atoms. Each actual Atom object has 
   * a pointer to its parent AtomArray and its array index in private 
   * members.  The public accessor method for each psuedo-member simply 
   * retrieves the rquired array element for this atom. In the current
   * implementation the position, force, atom type id, and isGhost 
   * flag are stored in true member variables of an Atom object, while
   * the velocity, mask, plan, and id (the global atom index) are all 
   * psuedo-members stored in separate arrays. See documentation of 
   * the private member localId_ and other comments in the Atom.h file 
   * for further implementation details.
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Atom
   {

   public:

      #ifdef UTIL_MPI
      /**
      * Return size of an atom packed into buffer for exchange, in bytes.
      */
      static int packedAtomSize();

      /**
      * Return size of ghost atom packed for communication, in bytes.
      */
      static int packedGhostSize();
      #endif

      /// \name Mutators
      //@{

      /**
      * Assignment.
      */
      Atom& operator = (const Atom& other);

      /**
      * Reset integer members to initial null values.
      */
      void clear();

      /**
      * Set unique global id for this Atom.
      *
      * This id must be unique among all atoms on all processors,
      * and must be less than totalAtomCapacity of the AtomStorage.
      *
      * \param Id unique atom id.
      */
      void setId(int Id);

      /**
      * Set the atom type index.
      *
      * \param Id integer index that identifies atom type
      */
      void setTypeId(int Id);

      /**
      * Mark as ghost or local atom for this processor.
      *
      * \param isGhost true if this is a ghost, or false if local.
      */
      void setIsGhost(bool isGhost);

      //@}
      /// \name Accessors (return non-const references)
      //@{

      /**
      * Get position Vector by reference.
      */
      Vector& position();

      /**
      * Get velocity Vector by reference.
      */
      Vector& velocity();

      /**
      * Get force Vector by reference.
      */
      Vector& force();

      /**
      * Get the associated Mask by reference.
      */
      Mask& mask();

      /**
      * Get communication plan by reference.
      */
      Plan& plan();

      #ifdef DDMD_MOLECULES
      /**
      * Get the AtomContext struct.
      *
      * A DdMd::AtomContext struct contains public members speciesId,
      * moleculeId and atomId that identify the species of molecule
      * to which this atom belongs, the index of the molecule within
      * it species, and the index of the atom with the molecule. 
      */
      AtomContext& context();
      #endif      

      //@}
      /// \name Accessors (return values and const references).
      //@{

      /**
      * Get unique global index for this atom.
      */
      int id() const;

      /**
      * Get atom type index.
      */
      int typeId() const;

      /**
      * Is this atom a ghost?
      */
      bool isGhost() const;

      /**
      * Get the position Vector (const reference).
      */
      const Vector& position() const;

      /**
      * Get the velocity Vector (const reference).
      */
      const Vector& velocity() const;

      /**
      * Get the force Vector (const reference).
      */
      const Vector& force() const;

      /**
      * Get the associated Mask by const reference.
      */
      const Mask& mask() const;

      /**
      * Get communication plan (const reference).
      */
      const Plan& plan() const;

      #ifdef DDMD_MOLECULES      
      /**
      * Get the context by reference (const reference).
      */
      const AtomContext& context() const;
      #endif

      //@}
      #if 0
      /// Get the shift IntVector by reference.
      IntVector& shift();

      /// Get the shift IntVector by const reference.
      const IntVector& shift() const;
      #endif

      #ifdef UTIL_MPI
      /// \name Pack and Unpack Methods (Interprocessor Communication)
      //@{

      /**
      * Pack an Atom into a send buffer, for exchange of ownership.
      *
      * Packs required data, increments buffer sendSize counter.
      *
      * \param buffer communication buffer
      */
      void packAtom(Buffer& buffer);

      /**
      * Unpack an atom from a recv buffer and receive ownership.
      *
      * \param buffer communication buffer
      */
      void unpackAtom(Buffer& buffer);

      /**
      * Pack a ghost Atom into a send buffer.
      *
      * Packs required data, increments buffer sendSize counter.
      *
      * \param buffer communication buffer
      */
      void packGhost(Buffer& buffer);

      /**
      * Unpack a ghost Atom from a recv buffer.
      *
      * Unpacks required data, decrements buffer recvSize counter.
      *
      * \param buffer communication buffer
      */
      void unpackGhost(Buffer& buffer);

      /**
      * Pack updated ghost position into send buffer.
      *
      * Packs position Vector, increments buffer sendSize counter.
      *
      * \param buffer communication buffer
      */
      void packUpdate(Buffer& buffer);

      /**
      * Unpack updated ghost position from recv buffer.
      *
      * Unpacks position Vector, decrements buffer recvSize counter.
      *
      * \param buffer communication buffer
      */
      void unpackUpdate(Buffer& buffer);

      /**
      * Pack update of ghost Atom force into send buffer.
      *
      * Packs force Vector, increments buffer sendSize counter.
      *
      * \param buffer communication buffer
      */
      void packForce(Buffer& buffer);

      /**
      * Unpack updated position of ghost Atom from recv buffer.
      *
      * Reads force from buffer and increments the atomic force
      * (rather than overwriting), then decrements recvSize counter.
      *
      * \param buffer communication buffer
      */
      void unpackForce(Buffer& buffer);

      //@}
      #endif

   private:

      /**
      * Position of atom.
      */
      Vector position_;

      /**
      * Integer index of atom type.
      */
      int typeId_;

      /**
      * Local id in Atom Array, set by AtomArray.
      *
      * The least signficant bit of this unsigned int stores a ghost flag,
      * which is 1 if this atom is a ghost, and 0 if it is a local atom. 
      * The remaining bits, which can be accessed as localId_ >> 1, store
      * the array index of this atom in the parent AtomArray. The array
      * index is set during allocation, and never changed thereafter.
      * The ghost flag may be changed at any time by the public function
      * void Atom::setIsGhost(bool).  
      */
      unsigned int localId_;

      /**
      * Force on atom.
      */
      Vector force_;

      /**
      * Pointer to parent AtomArray.
      */
      AtomArray* arrayPtr_;

      #ifdef UTIL_32BIT
      /**
      * On machines with 4 byte pointer, this pads the size to 64 bytes.
      */
      int pad_;
      #endif

      /**
      * Constructor.
      *
      * An Atom may be constructed only as an element of an AtomArray.
      * The AtomArray class is a friend of Atom, and may thus call this 
      * private constructor.
      */
      Atom();

      /**
      * Copy constructor.
      *
      * Private and not implemented to prohibit copy construction.
      */
      Atom(const Atom& other);

      friend class AtomArray;
      friend class Util::Memory;

   };

   // Inline methods

   /*
   * Set type Id for Atom.
   */
   inline void Atom::setTypeId(int typeId)
   {  typeId_ = typeId; }

   /*
   * Set isGhost flag.
   */
   inline void Atom::setIsGhost(bool isGhost)
   {
      if (isGhost) {
         // Set least significant bit of localId_ to 1
         localId_ = localId_ | 1;
      } else {
         // Set least significant bit of localId_ to 0
         localId_ = localId_ & ~1;
      }
   }

   /*
   * Get atom type Id.
   */
   inline int Atom::typeId() const
   {  return typeId_; }

   /*
   * Is this a ghost atom?
   */
   inline bool Atom::isGhost() const
   {
      // Return least significant bit of localId_
      return bool(localId_ & 1);
   }

   /*
   * Get position by reference.
   */
   inline Vector& Atom::position()
   {  return position_; }

   /*
   * Get position by const reference.
   */
   inline const Vector& Atom::position() const
   {  return position_; }

   /*
   * Get force by reference.
   */
   inline Vector& Atom::force()
   {  return force_; }

   /*
   * Get force by const reference.
   */
   inline const Vector& Atom::force() const
   {  return force_; }

   /*
   * Accessors for pseudo-members stored in separate arrays.
   *
   * An atom can be constructed only as an element of a parent AtomArray
   * (see documentation for private default constructor). Some data owned
   * by an Atom is stored in a set of private array that are allocated by
   * its AtomArray. In each such array, the element associated with this 
   * Atom is indexed by the localId_ of this atom, without the least 
   * signficant bit (see documentation of the localId_ member). The bit 
   * shift operation "localId_ >> 1" in each of the following functions 
   * strips off the least signficant bit, which is used to store the isGhost 
   * flag, and returns the relevant array index. Each arrays that stores
   * pseudo-member data is a private member of the parent AtomArray (e.g.,
   * AtomArray::velocities_, AtomArray::masks_, etc.) that is accessible
   * because Atom is a friend class of AtomArray.
   */

   /* 
   * Get reference to velocity.
   */
   inline Vector& Atom::velocity()
   {  return arrayPtr_->velocities_[localId_ >> 1]; }

   /*
   * Get velocity by const reference.
   */
   inline const Vector& Atom::velocity() const
   { return arrayPtr_->velocities_[localId_ >> 1]; }

   /*
   * Get the Mask by reference.
   */
   inline Mask& Atom::mask()
   {  return arrayPtr_->masks_[localId_ >> 1]; }

   /*
   * Get the Mask by const reference.
   */
   inline const Mask& Atom::mask() const
   {  return arrayPtr_->masks_[localId_ >> 1]; }

   /*
   * Get the communication plan by reference.
   */
   inline Plan& Atom::plan()
   {  return arrayPtr_->plans_[localId_ >> 1]; }

   /*
   * Get the communication plan by const reference.
   */
   inline const Plan& Atom::plan() const
   {  return arrayPtr_->plans_[localId_ >> 1]; }

   /*
   * Get the global id for this Atom.
   */
   inline int  Atom::id() const
   {  return arrayPtr_->ids_[localId_ >> 1]; }

   /*
   * Set the global id for this Atom.
   */
   inline void Atom::setId(int id)
   {  arrayPtr_->ids_[localId_ >> 1] = id; }

   #ifdef DDMD_MOLECULES   
   /* 
   * Get non-const reference to context.
   */
   inline AtomContext& Atom::context()
   {  return arrayPtr_->contexts_[localId_ >> 1]; }

   /*
   * Get context by const reference.
   */
   inline const AtomContext& Atom::context() const
   { return arrayPtr_->contexts_[localId_ >> 1]; }
   #endif
}
#endif
