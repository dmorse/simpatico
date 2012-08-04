#ifndef DDMD_MODIFIER_H
#define DDMD_MODIFIER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

#include <iostream>

namespace DdMd
{

   using namespace Util;
   class Simulation;

   /**
   * Abstract base for classes that can modify the integration loop.
   *
   * \ingroup DdMd_Modifier_Module
   */
   class Modifier : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Modifier(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Modifier();

      /**
      * Get interval value.
      */
      int interval() const;

      /**
      * Return true iff iStep is a multiple of the interval.
      *
      * \param iStep simulation step iStep
      */
      bool isAtInterval(long iStep) const;

      // Setup
   
      virtual void setupPostExchange(){};
      virtual void setupPostNeighbor(){};
      virtual void setupPostForce(){};
      virtual void setupEnd(){};
   
      // Integration

      virtual void preIntegrate() {};
      virtual void postIntegrate() {};
   
      virtual void preTransform() {};
      virtual void preExchange() {};
      virtual void postExchange() {};
      virtual void postNeighbor() {};
   
      virtual void preUpdate() {};
      virtual void postUpdate() {};
   
      virtual void preForce() {};
      virtual void postForce() {};
      virtual void endOfStep() {};
   
      // Communication
   
      virtual void pack_exchange(Buffer& buffer) {};
      virtual void unpack_exchange(Buffer& buffer) {};
      virtual void pack_update(Buffer& buffer) {};
      virtual void unpack_update(Buffer& buffer) {};
      virtual void pack_reverseUpdate(Buffer& buffer) {};
      virtual void unpack_reverseUpdate(Buffer& buffer) {};

   protected:

      /**
      * Read parameter interval from file.
      *
      * This function throws an exception if the value of interval
      * is not a multiple of Modifier::baseInterval, or if
      * baseInterval has not been set to a nonzero positive value.
      *
      * \param in input parameter file stream.
      */
      void readInterval(std::istream &in);

      /**
      * Get the parent Simulation by reference.
      */
      Simulation& simulation();

   private:

      /// Pointer to parent Simulation
      Simulation* simulationPtr_;

      /// Number of simulation steps between subsequent actions.
      long  interval_;

   };

   // Inline methods

   /*
   * Return interval value.
   */
   inline int Modifier::interval() const
   {  return interval_; }

   /*
   * Return true iff the iStep is a multiple of the interval.
   */
   inline bool Modifier::isAtInterval(long iStep) const
   {  return (iStep%interval_ == 0); }

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& Modifier::simulation()
   {  return *simulationPtr_; }

}
#endif
