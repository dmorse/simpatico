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

   /**
   * Abstract base for classes that can modify the integration loop.
   *
   * \ingroup DdMd_Action_Module
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

      // Setup
   
      virtual void setupPostExchange(){};
      virtual void setupPostNeighbor(){};
      virtual void setupPostForce(){};
   
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
   
      virtual void packExchange(Buffer& buffer) {};
      virtual void unpackExchange(Buffer& buffer) {};
      virtual void packUpdate(Buffer& buffer) {};
      virtual void unpackUpdate(Buffer& buffer) {};
      virtual void packreverseUpdate(Buffer& buffer) {};
      virtual void unpackreverseUpdate(Buffer& buffer) {};

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
      bool hasPackreverseUpdate();
      bool hasUnpackreverseUpdate();

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
      bool hasPackreverseUpdate_;
      bool hasUnpackreverseUpdate_;

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

   inline bool hasSetupPostExchange()
   {  return hasSetupPostExchange_; }

   inline bool hasSetupPostNeighbor()
   {  return hasSetupPostNeighbor_; }

   inline bool hasSetupPostForce()
   {  return hasSetupPostForce_; }

   inline bool hasPreIntegrate()
   {  return hasPreIntegrate_; }

   inline bool hasPostIntegrate()
   {  return hasPostIntegrate_; }

   inline bool hasPreTransform()
   {  return hasPreTransform_; }

   inline bool hasPreExchange()
   {  return hasPreExchange_; }

   inline bool hasPostExchange()
   {  return hasPostExchange_; }

   inline bool hasPostNeighbor()
   {  return hasPostNeighbor_; }

   inline bool hasPreUpdate()
   {  return hasPreUpdate_; }

   inline bool hasPostUpdate()
   {  return hasPostUpdate_; }

   inline bool hasPreForce()
   {  return hasPreForce_; }

   inline bool hasPostForce()
   {  return hasPostForce_; }

   inline bool hasEndOfStep()
   {  return hasEndOfStep_; }

   inline bool hasPackExchange()
   {  return hasUnpackExchange_; }

   inline bool hasPackUpdate()
   {  return hasUnpackUpdate_; }

   inline bool hasPackreverseUpdate()
   {  return hasUnpackreverseUpdate_; }

}
#endif
