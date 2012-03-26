#ifndef DDMD_ATOM_H
#define DDMD_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <ddMd/communicate/Plan.h>
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
   *
   * \ingroup DdMd_Chemistry_Module
   */
   class Atom
   {

   public:

      /**
      * Constructor.
      */
      Atom();
     
      // Use default destructor.
 
      /// \name Mutators
      //@{

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

      // Is this Atom a ghost? (0=false, 1=true)
      int isGhost_;

      /// Force on atom.
      Vector force_;                       

      /// Integer index for Atom within Simulation.
      int id_;                      

      // Mask (listed of bonded pairs)
      Mask mask_;   

      /// Atomic velocity.
      Vector velocity_;                       

      // Communication plan.
      Plan plan_;

      // IntVector shift_;  

   };

   // Inline methods

   // Constructor.
   inline Atom::Atom() :
     position_(0.0),
     typeId_(-1),
     isGhost_(0),
     force_(0.0),
     id_(-1),
     mask_(),
     velocity_(0.0),
     plan_()
   {}

   /**
   * Reset integer members to null values.
   */
   inline void Atom::clear()
   {
      typeId_ = -1;
      id_ = -1;
      isGhost_ = 0;
      mask_.clear();
      plan_.setFlags(0);
   }

   // Set unique global index for Atom.
   inline void Atom::setId(int id) 
   {  id_ = id; }

   // Set type Id for Atom.
   inline void Atom::setTypeId(int typeId) 
   {  typeId_ = typeId; }

   // Set type Id for Atom.
   inline void Atom::setIsGhost(bool isGhost) 
   {  isGhost_ = isGhost ? 1 : 0; }

   // Get global id for Atom.
   inline int  Atom::id() const
   {  return id_; }

   // Get type Id.
   inline int Atom::typeId() const
   {  return typeId_; }

   // Is this a ghost atom?
   inline bool Atom::isGhost() const
   {  return bool(isGhost_); }

   // Get reference to communication plan.
   inline Plan& Atom::plan()
   {  return plan_; }

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

   // Get const reference to communication plan.
   inline const Plan& Atom::plan() const
   {  return plan_; }

   // Get the associated mask.
   inline Mask& Atom::mask() 
   {  return mask_; }

   // Get a const reference to the mask.
   inline const Mask& Atom::mask() const
   {  return mask_; }

}
#endif
