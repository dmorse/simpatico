#ifndef ATOM_H
#define ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#include "Mask.h"
#include <util/space/Vector.h>

namespace DdMd
{

   using namespace Util;

   /**
   * A point particle.
   *
   * Each Atom has position, vector, and force Vector objects, an integer 
   * atom type Id, and a global integer id.
   * Each Atom has an associated Mask object. 
   */
   class Atom
   {

   public:

      // Default constructor.
      Atom();
      
      /// \name Mutators
      //@{

      /**
      * Set unique id.
      *  
      * \param Id unique atom id.
      */
      void setId(int Id);

      /**
      * Set the atomic type index.
      *  
      * \param Id integer index that identifies atom type
      */
      void setTypeId(int Id);

      /**
      * Set rank of the processor that owns this Atom.
      *  
      * \param ownerRank rank of parent processor.
      */
      void setOwnerRank(int ownerRank);

      /**
      * Set the atomic type index.
      *  
      * \param isGhost true if this is a ghost, or false if local.
      */
      void setIsGhost(bool isGhost);

      /**
      * Set or unset send marker.
      *  
      * \param sendMark true to mark for sending, false otherwise.
      */
      void setSendMark(bool sendMark);

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

      //@}
      /// \name Accessors 
      //@{

      /// Get unique global index.
      int   id() const;

      /// Get atom type index.
      int   typeId() const;

      /// Get rank of processor that owns this atom.
      int   ownerRank() const;

      /// Is this atom a ghost?
      bool  isGhost() const;

      /// Is this atom marked for sending?
      bool  sendMark() const;

      /// Get the position Vector (const reference).
      const Vector& position() const;

      /// Get the velocity Vector (const reference).
      const Vector& velocity() const;

      /// Get the force Vector (const reference).
      const Vector& force() const;

      //@}

      #if 0
      /// Get the associated Mask by reference.
      Mask& mask();

      /// Get the associated Mask by const reference.
      const Mask& mask() const;

      /// Get the shift IntVector by reference.
      IntVector& shift();

      /// Get the shift IntVector by const reference.
      const IntVector& shift() const;
      #endif

   private:

      /// Position of atom.
      Vector position_;                       

      /// Integer index of atom type.
      int  typeId_;                         

      /// Integer index for Atom within Simulation.
      int  id_;                      

      /// Force on atom.
      Vector force_;                       

      // Is this Atom a ghost (0=false, 1=true)
      int  isGhost_;

      // Is this Atom a ghost (0=false, 1=true)
      int  ownerRank_;

      // Integer index of molecule
      // int  moleculeId_;                      

      /// Atomic velocity.
      Vector velocity_;                       

      // Is this Atom a ghost (0=false, 1=true)
      int  sendMark_;

      // Mask      mask_;   
      // IntVector shift_;  

   };

   // Inline methods

   // Constructor.
   inline Atom::Atom() :
     position_(0.0),
     typeId_(-1),
     id_(-1),
     force_(0.0),
     isGhost_(0),
     ownerRank_(-1),
     // moleculeId(-1),
     velocity_(0.0),
     sendMark_(0)
   {}

   // Set unique global index for Atom.
   inline void Atom::setId(int id) 
   {  id_ = id; }

   // Set type Id for Atom.
   inline void Atom::setTypeId(int typeId) 
   {  typeId_ = typeId; }

   // Set rank of owner processor.
   inline void Atom::setOwnerRank(int ownerRank) 
   {  ownerRank_ = ownerRank; }

   // Set type Id for Atom.
   inline void Atom::setIsGhost(bool isGhost) 
   {  isGhost_ = isGhost ? 1 : 0; }

   // Set type Id for Atom.
   inline void Atom::setSendMark(bool sendMark) 
   {  sendMark_ = sendMark ? 1 : 0; }

   // Get global id for Atom.
   inline int  Atom::id() const
   {  return id_; }

   // Get type Id.
   inline int Atom::typeId() const
   {  return typeId_; }

   // Get type Id.
   inline int Atom::ownerRank() const
   {  return ownerRank_; }

   // Is this a ghost atom?
   inline bool Atom::isGhost() const
   {  return bool(isGhost_); }

   // Is this atom marked for sending?
   inline bool Atom::sendMark() const
   {  return bool(sendMark_); }

   // Get reference to position.
   inline Vector& Atom::position()
   {  return position_; }

   // Get reference to velocity.
   inline Vector& Atom::velocity()
   {  return velocity_; }

   // Get reference to force.
   inline Vector& Atom::force()
   {  return force_; }

   // Get const reference to position.
   inline const Vector& Atom::position() const
   {  return position_; }

   // Get const reference to velocity.
   inline const Vector& Atom::velocity() const
   {  return velocity_; }

   // Get const reference to force.
   inline const Vector& Atom::force() const
   {  return force_; }

   #if 0
   // Get the associated mask.
   inline Mask& Atom::mask() 
   {  return masks_[id_]; }

   // Get a const reference to the mask.
   inline const Mask& Atom::mask() const
   {  return masks_[id_]; }
   #endif

}
#endif
