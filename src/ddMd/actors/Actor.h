#ifndef DDMD_ACTOR_H
#define DDMD_ACTOR_H

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/Bit.h>              // constants in namespace

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace DdMd
{

   using namespace Util;
   class Simulation;

   /**
   * An Actor defines actions taken during integration.
   *
   * The Actor base class declares a variety of virtual "action" methods 
   * that may be re-implemented by subclasses to define actions that 
   * should be executed at specific points before, during, or after the 
   * main integration loop. The name of each such action method describes 
   * when it will be invoked, if activated (see below). For example, the 
   * postForce() method is invoked within the main integration loop 
   * immediately after evaluation of all forces. 
   *
   * For each action method, there is a boolean flag. Each action method 
   * will be executed if and only if the corresponding flag is set to 
   * true. All flags are set to false in the base class constructor, but
   * should be set to true in the subclass constructor to activate a
   * method. 
   *
   * Each Actor class also has an associated interval integer member. 
   * An activated action method of an Actor will be executed only on 
   * time steps that are multiples of this interval. Different intervals 
   * are not defined for different action methods: a single interval 
   * value is defined for each Actor object. The value of the interval 
   * is initialized to 1 (every time step) in the base Actor class 
   * constructor, but may be reset to a greater value in the subclass 
   * readParam method, by calling the protected readInterval() method.
   */
   class Actor : public ParamComposite
   {

   public:

      /**
      * Default constructor (for unit testing)
      */
      Actor();

      /**
      * Constructor (for use in simulation).
      */
      Actor(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Actor();

      /// \name Setup actions 
      //@{ 
   
      virtual void setupPostExchange(){};
      virtual void setupPostNeighbor(){};
      virtual void setupPostForce(){};
  
      //@} 
      /// \name Integration actions (within main loop)
      //@{ 

      virtual void preIntegrate1() {};
      virtual void postIntegrate1() {};
   
      virtual void preTransform() {};
      virtual void preExchange() {};
      virtual void postExchange() {};
      virtual void postNeighbor() {};
      virtual void preUpdate() {};
      virtual void postUpdate() {};
      virtual void preForce() {};
      virtual void postForce() {};
      virtual void endOfStep() {};

      //@} 
      /// \name Interprocessor communication actions 
      //@{ 

      virtual void packExchange() {};
      virtual void unpackExchange() {};
      virtual void packUpdate() {};
      virtual void unpackUpdate() {};
      virtual void packReverseUpdate() {};
      virtual void unpackReverseUpdate() {};

      //@} 
      /// \name Bit Flags 
      //@{

      /**
      * Bit flag constant associated with particular virtual methods.
      */
      class Flags 
      { 
      public:
         static const Bit SetupPostExchange;
         static const Bit SetupPostNeighbor;
         static const Bit SetupPostForce;
         static const Bit PreIntegrate1;
         static const Bit PostIntegrate1;
         static const Bit PreTransform;
         static const Bit PreExchange;
         static const Bit PostExchange;
         static const Bit PostNeighbor;
         static const Bit PreUpdate;
         static const Bit PostUpdate;
         static const Bit PreForce;
         static const Bit PostForce;
         static const Bit EndOfStep;
         static const Bit Exchange;
         static const Bit Update;
         static const Bit ReverseUpdate;
      };

      /**
      * Return true if a flag is set, false otherwise.
      */
      bool isSet(Bit flag);

      //@} 
      /// \name Interval methods
      //@{

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

      //@} 

   protected:

      /**
      * Read parameter interval from file.
      *
      * This function throws an exception if the value of interval
      * is not a multiple of Actor::baseInterval, or if
      * baseInterval has not been set to a nonzero positive value.
      *
      * \param in input parameter file stream.
      */
      void readInterval(std::istream &in);

      /**
      * Get the parent Simulation by reference.
      */
      Simulation& simulation();

      /**
      * Set the flag associated with a particular action.
      */
      void set(Bit flag);

   private:

      unsigned int flags_;

      /// Number of simulation steps between subsequent actions.
      long  interval_;

      /// Pointer to parent Simulation
      Simulation* simulationPtr_;

   };

   // Inline methods

   /*
   * Return interval value.
   */
   inline int Actor::interval() const
   {  return interval_; }

   /*
   * Return true iff the iStep is a multiple of the interval.
   */
   inline bool Actor::isAtInterval(long iStep) const
   {  return (iStep%interval_ == 0); }

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& Actor::simulation()
   {  return *simulationPtr_; }

}
#endif
