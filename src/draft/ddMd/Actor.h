#ifndef DDMD_ACTOR_H
#define DDMD_ACTOR_H

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
   * An actor defines actions taken during execution.
   *
   * Methods of an Actor implement actions that are executed at
   * specified points during execution of a simulation. The Actor 
   * base class declares a variety of virtual "action" methods that 
   * may be re-implemented by subclasses to define actions that should 
   * be executed at specific points before, during, or after the main 
   * integration loop. The name of each such action method describes 
   * when it will be invoked. For example, the postForce() method is 
   * invoked within the main integration loop immediately after 
   * evaluation of all forces. 
   *
   * For each action method, there is a associated boolean data member 
   * whose name is given by the name of the associated method prefixed 
   * by "has", and suffixed by an underscore. For example, the boolean
   * flag associated with the postForce() method is hasPostForce_.
   * Each such boolean member is initialized to false in the constructor
   * of the Actor base class, and must be reset to true by the constructor
   * of a subclass to activate execution of the associated method. Each
   * action method will be executed if and only if the corresponding 
   * boolean flag is set to true in a subclass constructor. Each action
   * method has an empty default implementation, so activating a method
   * without re-implementing the method would do nothing.
   *
   * Each subclass of Actor must thus re-implement one or more of the
   * virtual action methods, and set the associated boolean flags to 
   * true to activate these methods. For example, to implement only
   * one action that should be executed immediately after forces are 
   * evaluated, one would create a a subcalss of Actor that:
   *
   * (1) Re-implements the virtual postForce() method so as to define
   * the desired action, and 
   *
   * (2) Re-initializes the protected variable hasPostForce_ to true 
   * in the constructor for that subclass. 
   *
   * Each Actor class also has an associated interval integer member. 
   * A re-implemented action method of an Actor will be executed only 
   * on time steps that are multiples of this interval. Different 
   * intervals are not defined for different action methods: a single
   * interval value is defined for each Actor object. The value of the 
   * interval is initialized to 1 (every time step) in the base Actor 
   * class constructor, but may be reset to a greater value in the 
   * subclass readParam method, by calling the readInterval() method.
   */
   class Actor : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Actor(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Actor();

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

      /// \name Pre-integration actions 
      //@{ 
   
      virtual void setupPostExchange(){};
      virtual void setupPostNeighbor(){};
      virtual void setupPostForce(){};
  
      //@} 
      /// \name Integration actions (within main loop)
      //@{ 

      virtual void preIntegrate() {};
      virtual void postIntegrate() {};
   
      virtual void preTransform() {};
      virtual void preExchange() {};
      virtual void packExchange() {};
      virtual void unpackExchange() {};
      virtual void postExchange() {};
      virtual void postNeighbor() {};
   
      virtual void preUpdate() {};
      virtual void packUpdate() {};
      virtual void unpackUpdate() {};
      virtual void postUpdate() {};
   
      virtual void preForce() {};
      virtual void packReverseUpdate() {};
      virtual void unpackReverseUpdate() {};
      virtual void postForce() {};

      virtual void endOfStep() {};
   

      //@} 
      /// \name Boolean flags (accessors)
      //@{ 

      bool hasSetupPostExchange();
      bool hasSetupPostNeighbor();
      bool hasSetupPostForce();

      bool hasPreIntegrate();
      bool hasPostIntegrate();
      bool hasPreTransform();
      bool hasPreExchange();
      bool hasPostExchange();
      bool hasPostNeighbor();
      bool hasPreUpdate();
      bool hasPostUpdate();
      bool hasPreForce();
      bool hasPostForce();
      bool hasEndOfStep();

      bool hasPackExchange();
      bool hasUnpackExchange();
      bool hasPackUpdate();
      bool hasUnpackUpdate();
      bool hasPackReverseUpdate();
      bool hasUnpackReverseUpdate();

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

      /// Number of simulation steps between subsequent actions.
      long  interval_;

      /// \name Boolean flags 
      //@{ 

      bool hasSetupPostExchange_;
      bool hasSetupPostNeighbor_;
      bool hasSetupPostForce_;
      bool hasPreIntegrate_;
      bool hasPostIntegrate_;
      bool hasPreTransform_;
      bool hasPreExchange_;
      bool hasPostExchange_;
      bool hasPostNeighbor_;
      bool hasPreUpdate_;
      bool hasPostUpdate_;
      bool hasPreForce_;
      bool hasPostForce_;
      bool hasEndOfStep_;
      bool hasPackExchange_;
      bool hasUnpackExchange_;
      bool hasPackUpdate_;
      bool hasUnpackUpdate_;
      bool hasPackReverseUpdate_;
      bool hasUnpackReverseUpdate_;

      //@}

   private:

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

   inline bool Actor::hasSetupPostExchange()
   {  return hasSetupPostExchange_; }

   inline bool Actor::hasSetupPostNeighbor()
   {  return hasSetupPostNeighbor_; }

   inline bool Actor::hasSetupPostForce()
   {  return hasSetupPostForce_; }

   inline bool Actor::hasPreIntegrate()
   {  return hasPreIntegrate_; }

   inline bool Actor::hasPostIntegrate()
   {  return hasPostIntegrate_; }

   inline bool Actor::hasPreTransform()
   {  return hasPreTransform_; }

   inline bool Actor::hasPreExchange()
   {  return hasPreExchange_; }

   inline bool Actor::hasPostExchange()
   {  return hasPostExchange_; }

   inline bool Actor::hasPostNeighbor()
   {  return hasPostNeighbor_; }

   inline bool Actor::hasPreUpdate()
   {  return hasPreUpdate_; }

   inline bool Actor::hasPostUpdate()
   {  return hasPostUpdate_; }

   inline bool Actor::hasPreForce()
   {  return hasPreForce_; }

   inline bool Actor::hasPostForce()
   {  return hasPostForce_; }

   inline bool Actor::hasEndOfStep()
   {  return hasEndOfStep_; }

   inline bool Actor::hasPackExchange()
   {  return hasPackExchange_; }

   inline bool Actor::hasUnpackExchange()
   {  return hasUnpackExchange_; }

   inline bool Actor::hasPackUpdate()
   {  return hasPackUpdate_; }

   inline bool Actor::hasUnpackUpdate()
   {  return hasUnpackUpdate_; }

   inline bool Actor::hasPackReverseUpdate()
   {  return hasPackReverseUpdate_; }

   inline bool Actor::hasUnpackReverseUpdate()
   {  return hasUnpackReverseUpdate_; }

}
#endif
