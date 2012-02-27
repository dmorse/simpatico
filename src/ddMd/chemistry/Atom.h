#ifndef ATOM_H
#define ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/communicate/Plan.h>
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
      * Set the atomic type index.
      *  
      * \param isGhost true if this is a ghost, or false if local.
      */
      void setIsGhost(bool isGhost);

      /**
      * Set the postMark.
      *  
      * \param postMark true if this is marked for sending, false otherwise.
      */
      void setPostMark(bool postMark);

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

      #if 0
      /**
      * Get communication plan by reference.
      */
      Plan& plan();
      #endif

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

      #if 0
      /// Communication plan (const reference).
      const Plan& plan() const;
      #endif

      /// Return postMark (true if marked for sending).
      bool postMark() const;

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
      int typeId_;                         

      /// Integer index for Atom within Simulation.
      int id_;                      

      /// Force on atom.
      Vector force_;                       

      // Is this Atom a ghost? (0=false, 1=true)
      int isGhost_;

      // Is this Atom marked for sending? (0=false, 1=true)
      int postMark_;

      #if 0
      // Is this Atom a ghost (0=false, 1=true)
      Plan plan_;
      #endif

      /// Atomic velocity.
      Vector velocity_;                       

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
     postMark_(0),
     //plan_(),
     velocity_(0.0)
   {}

   // Set unique global index for Atom.
   inline void Atom::setId(int id) 
   {  id_ = id; }

   // Set type Id for Atom.
   inline void Atom::setTypeId(int typeId) 
   {  typeId_ = typeId; }

   // Set type Id for Atom.
   inline void Atom::setIsGhost(bool isGhost) 
   {  isGhost_ = isGhost ? 1 : 0; }

   // Set type Id for Atom.
   inline void Atom::setPostMark(bool postMark) 
   {  postMark_ = postMark ? 1 : 0; }

   // Get global id for Atom.
   inline int  Atom::id() const
   {  return id_; }

   // Get type Id.
   inline int Atom::typeId() const
   {  return typeId_; }

   // Is this a ghost atom?
   inline bool Atom::isGhost() const
   {  return bool(isGhost_); }

   // Is this atom marked for sendng?
   inline bool Atom::postMark() const
   {  return bool(postMark_); }

   #if 0
   // Get reference to communication plan.
   inline Plan& Atom::plan()
   {  return plan_; }
   #endif

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
   // Get const reference to communication plan.
   inline const Plan& Atom::plan() const
   {  return plan_; }
   #endif

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
