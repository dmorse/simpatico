#ifndef DDMD_ATOM_EVENTS_H
#define DDMD_ATOM_EVENTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Veronica Chappa and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   /**
   * Base class for events signalling change in atom data.
   */
   class AtomEvent
   {

   public:

      unsigned int flags()
      {  return flags_; }

   protected:

      /// Constructor (constructed to prevent direct instantiation).
      AtomAddEvent(flags)
       : flags_(flags)
      {}

   private:

      /// Integer bitfield of flags.
      unsigned int flags_;

   };

   /**
   * Event signalling change in atom positions.
   */
   class PositionEvent
   {
   public:

      /// Constructor.
      PositionEvent()
       : AtomEvent(1)
      {}

   };

   /**
   * Event signalling change in atom velocities.
   */
   class VelocityEvent
   {
   public:

      /// Constructor.
      VelocityEvent()
       : AtomEvent(2)
      {}

   };

   /**
   * Event signalling change in atom forces.
   */
   class ForceEvent
   {
   public:

      /// Constructor.
      ForceEvent()
       : AtomEvent(4)
      {}

   };

} 
#endif
