#ifndef DDMD_ATOM_H
#define DDMD_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#define UTIL_32BIT

#include <util/space/Vector.h>
#include <ddMd/chemistry/Mask.h>
#include <ddMd/communicate/Plan.h>
#include "AtomArray.h"

namespace DdMd
{

   using namespace Util;

   /**
   * A point particle.
   *
   * Each Atom has:
   *
   *   - a position Vector
   *   - a force  Vector
   *   - a velocity Vector
   *   - an integer atom type Id
   *   - a boolean isGhost flag
   *   - a global integer id
   *   - a Mask (list of other atoms with masked pair interactions)
   *   - a communication plan
   *
   * An Atom may only be constructed as an element of an AtomArray.
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Atom
   {

   public:

      /// \name Mutators
      //@{

      /**
      * Assignment.
      */
      Atom& operator= (const Atom& other);

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

      //@}
      /// \name Accessors (return values and const references).
      //@{

      /// Get unique global atom index.
      int  id() const;

      /// Get atom type index.
      int  typeId() const;

      /// Is this atom a ghost?
      bool isGhost() const;

      /// Get the position Vector (const reference).
      const Vector& position() const;

      /// Get the velocity Vector (const reference).
      const Vector& velocity() const;

      /// Get the force Vector (const reference).
      const Vector& force() const;

      /// Get communication plan (const reference).
      const Plan& plan() const;

      /// Get the associated Mask by const reference.
      const Mask& mask() const;

      //@}

      #if 0
      /// Get the shift IntVector by reference.
      IntVector& shift();

      /// Get the shift IntVector by const reference.
      const IntVector& shift() const;
      #endif

   private:

      /*
      * Position of atom.
      */
      Vector position_;

      /*
      * Integer index of atom type.
      */
      int typeId_;

      /*
      * Local id in Atom Array, set by AtomArray.
      * Least signicant bit: Is this Atom a ghost? (0=false, 1=true).
      * Only the ghost bit may be changed after allocation of AtomArray
      */
      unsigned int localId_;

      /*
      * Force on atom.
      */
      Vector force_;

      /*
      * Pointer to parent AtomArray.
      */
      AtomArray* arrayPtr_;

      #ifdef UTIL_32BIT
      /*
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
   * Accessors for associated data stored in arrays allocated  by AtomArray.
   *
   * An atom can be constructed only as an element of a parent AtomArray
   * (see documentation for private default constructor). Some data owned
   * by an Atom is stored in a set of private array that are allocated by
   * its AtomArray. In each such array, the element associated with this 
   * Atom is indexed by the localId_ of this atom, without the isGhost
   * flag. The bit shift operation "localId_ >> 1" in each of the following
   * functions strips off the least signficant bit, which is used to store
   * the isGhost flag, and returns the relevant array index. This array index
   * is set for each atom when AtomArray is allocated, and is not changed
   * thereafter.
   */
   
   // Get reference to velocity.
   inline Vector& Atom::velocity()
   {  return arrayPtr_->velocities_[localId_ >> 1]; }

   // Get velocity by const reference.
   inline const Vector& Atom::velocity() const
   { return arrayPtr_->velocities_[localId_ >> 1]; }

   // Get the Mask by reference.
   inline Mask& Atom::mask()
   {  return arrayPtr_->masks_[localId_ >> 1]; }

   // Get the Mask by const reference.
   inline const Mask& Atom::mask() const
   {  return arrayPtr_->masks_[localId_ >> 1]; }

   // Get the communication plan by reference.
   inline Plan& Atom::plan()
   {  return arrayPtr_->plans_[localId_ >> 1]; }

   // Get the communication plan by const reference.
   inline const Plan& Atom::plan() const
   {  return arrayPtr_->plans_[localId_ >> 1]; }

   // Get the global id for this Atom.
   inline int  Atom::id() const
   {  return arrayPtr_->ids_[localId_ >> 1]; }

   // Set the global id for this Atom.
   inline void Atom::setId(int id)
   {  arrayPtr_->ids_[localId_ >> 1] = id; }

}
#endif
