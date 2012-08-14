#ifndef DDMD_ACTION_H
#define DDMD_ACTION_H

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
   class Buffer;

   /**
   * ABC for actins during integration loop.
   *
   * Subclasses of Action implement actions that should be executed
   * at specified points in the main integration loop. Each subclass
   * of action should overload one or more of the virtual methods 
   * with empty default implementations to execute actions during 
   * setup, integration, or communication, for which the name of 
   * the method indicates where the corresponding action should be 
   * executed.  For each such re-implemented method, the constructor 
   * of the subclass must also set the corresponding protected bool 
   * variable to true.
   *
   * For example, to implement an action should be executed 
   * immediately after all forces are evaluated, a developer must
   * create a subcalss of Action that:
   *
   * (1) Re-implements the virtual postForce() method so as to 
   * execute the desired action, and 
   *
   * (2) Reinitializes the protected variable hasPostForce_ to 
   * true in the constructor for that subclass. 
   * 
   * Each of the re-implemented methods of an Action are executed
   * only if the corresponding bool variable is also set to true,
   * and only on time steps that are multiples of an associated
   * interval. A single interval is defined for the entire class, 
   * rather than for each re-implemented method. The value of the
   * interval should set by the readParam method, by calling the
   * protected readInterval method. 
   */
   class Action : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Action(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Action();

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

      // Setup (called before main loop)
   
      virtual void setupPostExchange(){};
      virtual void setupPostNeighbor(){};
      virtual void setupPostForce(){};
   
      // Integration (called within the main loop)

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
   
      // Communication (called within methods of Exchanger)
   
      virtual void packExchange(Buffer& buffer) {};
      virtual void unpackExchange(Buffer& buffer) {};
      virtual void packUpdate(Buffer& buffer) {};
      virtual void unpackUpdate(Buffer& buffer) {};
      virtual void packReverseUpdate(Buffer& buffer) {};
      virtual void unpackReverseUpdate(Buffer& buffer) {};

      // Boolean accessors

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

   protected:

      /**
      * Read parameter interval from file.
      *
      * This function throws an exception if the value of interval
      * is not a multiple of Action::baseInterval, or if
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

      // Boolean flag members

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

   private:

      /// Pointer to parent Simulation
      Simulation* simulationPtr_;

   };

   // Inline methods

   /*
   * Return interval value.
   */
   inline int Action::interval() const
   {  return interval_; }

   /*
   * Return true iff the iStep is a multiple of the interval.
   */
   inline bool Action::isAtInterval(long iStep) const
   {  return (iStep%interval_ == 0); }

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& Action::simulation()
   {  return *simulationPtr_; }

   inline bool Action::hasSetupPostExchange()
   {  return hasSetupPostExchange_; }

   inline bool Action::hasSetupPostNeighbor()
   {  return hasSetupPostNeighbor_; }

   inline bool Action::hasSetupPostForce()
   {  return hasSetupPostForce_; }

   inline bool Action::hasPreIntegrate()
   {  return hasPreIntegrate_; }

   inline bool Action::hasPostIntegrate()
   {  return hasPostIntegrate_; }

   inline bool Action::hasPreTransform()
   {  return hasPreTransform_; }

   inline bool Action::hasPreExchange()
   {  return hasPreExchange_; }

   inline bool Action::hasPostExchange()
   {  return hasPostExchange_; }

   inline bool Action::hasPostNeighbor()
   {  return hasPostNeighbor_; }

   inline bool Action::hasPreUpdate()
   {  return hasPreUpdate_; }

   inline bool Action::hasPostUpdate()
   {  return hasPostUpdate_; }

   inline bool Action::hasPreForce()
   {  return hasPreForce_; }

   inline bool Action::hasPostForce()
   {  return hasPostForce_; }

   inline bool Action::hasEndOfStep()
   {  return hasEndOfStep_; }

   inline bool Action::hasPackExchange()
   {  return hasPackExchange_; }

   inline bool Action::hasUnpackExchange()
   {  return hasUnpackExchange_; }

   inline bool Action::hasPackUpdate()
   {  return hasPackUpdate_; }

   inline bool Action::hasUnpackUpdate()
   {  return hasUnpackUpdate_; }

   inline bool Action::hasPackReverseUpdate()
   {  return hasPackReverseUpdate_; }

   inline bool Action::hasUnpackReverseUpdate()
   {  return hasUnpackReverseUpdate_; }

}
#endif
