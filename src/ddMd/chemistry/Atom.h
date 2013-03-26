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
   * Each Atom has position, vector, and force Vector objects, an integer
   * atom type Id, and a global integer id.
   * Each Atom has an associated Mask object.
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
      * Set unique id for this Atom.
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
      * Mark as ghost or local atom.
      *
      * \param isGhost true if this is a ghost, or false if local.
      */
      void setIsGhost(bool isGhost);

      //@}
      /// \name Accessors (non-const references)
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

      /// Get the associated Mask by reference.
      Mask& mask();

      /**
      * Get communication plan by reference.
      */
      Plan& plan();

      //@}
      /// \name Accessors
      //@{

      /// Get unique global index.
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

      /// Position of atom.
      Vector position_;

      /// Integer index of atom type.
      int typeId_;

      /// Least signicant bit: Is this Atom a ghost? (0=false, 1=true)
      /// Remaining bits: Local id in Atom Array, set by AtomArray.
      unsigned int localId_;

      /// Force on atom.
      Vector force_;

      /// Pointer to atom array
      AtomArray* arrayPtr_;

      #ifdef UTIL_32BIT
      int pad_;
      #endif

      // IntVector shift_;

      /**
      * Constructor (called by AtomArray).
      */
      Atom();

      /**
      * Copy constructor (not implemented).
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
         // Set least significant bit to 1
         localId_ = localId_ | 1;
      } else {
         // Set least significant bit to 0
         localId_ = localId_ & ~1;
      }
   }

   // Get type Id.
   inline int Atom::typeId() const
   {  return typeId_; }

   // Is this a ghost atom?
   inline bool Atom::isGhost() const
   {
      // Return least significant bit
      return bool(localId_ & 1);
   }

   // Get reference to position.
   inline Vector& Atom::position()
   {  return position_; }

   // Get const reference to position.
   inline const Vector& Atom::position() const
   {  return position_; }

   // Get reference to force.
   inline Vector& Atom::force()
   {  return force_; }

   // Get const reference to force.
   inline const Vector& Atom::force() const
   {  return force_; }

   // Indirect access via other arrays
   
   // Get const reference to velocity.
   inline const Vector& Atom::velocity() const
   {
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->velocities_[localId_ >> 1]; 
   }

   // Get reference to velocity.
   inline Vector& Atom::velocity()
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->velocities_[localId_ >> 1]; 
   }

   // Get the associated mask.
   inline Mask& Atom::mask()
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->masks_[localId_ >> 1]; 
   }

   // Get a const reference to the mask.
   inline const Mask& Atom::mask() const
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->masks_[localId_ >> 1]; 
   }

   // Get reference to communication plan.
   inline Plan& Atom::plan()
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->plans_[localId_ >> 1]; 
   }

   // Get const reference to communication plan.
   inline const Plan& Atom::plan() const
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->plans_[localId_ >> 1]; 
   }

   // Get global id for Atom.
   inline int  Atom::id() const
   {  
      //int j = (int)(localId_ >> 1);  
      return arrayPtr_->ids_[localId_ >> 1]; 
   }

   // Set unique global index for Atom.
   inline void Atom::setId(int id)
   {  
      //int j = (int)(localId_ >> 1);  
      arrayPtr_->ids_[localId_ >> 1] = id; 
   }

}
#endif
