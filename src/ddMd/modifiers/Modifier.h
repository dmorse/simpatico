#ifndef DDMD_MODIFIER_H
#define DDMD_MODIFIER_H

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

   class Simulation;
   using namespace Util;

   /**
   * A Modifier can modify the Simulation during MD integration.
   *
   * The Modifier base class declares a variety of virtual "action" methods 
   * that may be re-implemented by subclasses to define actions that can
   * modify the system, and that, if activated, should be executed as specific 
   * points during, or after the main integration loop. The name of each such 
   * action method describes when it will be invoked if activated. For example, 
   * the postForce() method is invoked within the main integration loop 
   * immediately after evaluation of all forces. 
   *
   * Each action method may be activated by setting a corresponding flag.
   * Each action method is executed if and only if the corresponding flag
   * has is set.  All flags are cleared in the base class constructor.
   * Subclasses that implement specific methods must set the corresponding
   * flag for each such method in the subclass constructor to activate all 
   * re-implemented methods. 
   *
   * Each Modifier also has an associated interval integer member. 
   * An activated action method of an Modifier will be executed only on 
   * time steps that are multiples of this interval. Different intervals 
   * are not defined for different action methods: a single interval 
   * value is defined for each Modifier object. The value of the interval 
   * is initialized to 1 (every time step) in the base Modifier class 
   * constructor, but may be reset to a greater value in the subclass 
   * readParam method, by calling the protected readInterval() method.
   *
   * The design of the Modifer class is inspired by the Lammps "Fix" 
   * class, which it closely resembles. 
   *
   * \ingroup DdMd_Modifier_Module
   */
   class Modifier : public ParamComposite
   {

   public:

      /**
      * Default constructor (for unit testing)
      */
      Modifier();

      /**
      * Constructor (for use in simulation).
      */
      Modifier(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Modifier();

      /// \name Setup and integration actions 
      //@{ 
   
      virtual void setup(){};
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

      #if 1
      /**
      * Bit flag constants associated with particular actions.
      *
      * The flag for each reimplemented method should be set 
      * in the subclass constructor by passing the appropriate
      * Bit constant to the protected void Modifier::set(Bit) 
      * function. 
      *
      * For example, to activate an action that should be
      * invoked at the end of each time step, one would call
      * \code
      *    set(Flags::EndOfStep);
      * \endcode
      * within the body of the constructor for the subclass
      * of Modifer that defines this member action function. 
      */
      class Flags 
      { 
      public:
         static const Bit Setup;
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
      #endif

      /**
      * Return true if a flag is set, false otherwise.
      */
      bool isSet(Bit flag) const;

      /**
      * Return unsigned int representation of all bit flags.
      */
      unsigned int flags() const;

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
