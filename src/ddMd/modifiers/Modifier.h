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
   * modify the system, each of which will be executed at a specific points 
   * during or after the main integration loop. The name of each such action
   * method describes when it will be invoked. For example, the postForce() 
   * method will be invoked (if at all) immediately after evaluation of all 
   * forces within the main integration loop. The Modifier base class 
   * buffer declares several virtual communication methods that may be 
   * reimplemented by sublcasses to pack additional information into the
   * messages that are communicated between processors during exchange and
   * update steps.
   *
   * Each subclass of Modifier will generally implement only a few of many
   * possible action and communication methods. To record which such methods 
   * have been implemented, each subclass must "activate" each re-implemented
   * method by invoking the set() method within its constructor, and 
   * passing set() a named constant that identifies the name of the method.
   * Thus, for example, a sublclass that provides an implementation of the
   * virtual postForce() function should also invoke set(Flags::PostForce) 
   * in its constructor. The constant Flags::PostForce and other named 
   * constants are defined as static constant members of the Modifier::Flags 
   * nested class. The name of the constant that is used to activate a
   * corresponding virtual method is generally a capitalized version of 
   * the name of the corresponding method.
   *
   * Each Modifier also has an associated interval integer member. An
   * activated action method of a Modifier will be executed only on time 
   * steps that are multiples of this interval. Each Modifier object has 
   * a single interval value that applies to all action methods that it
   * implements.  The value of the interval is initialized to 1 (once per 
   * time step) in the base Modifier class constructor, but may be reset 
   * in the subclass readParameters method by invoking the readInterval() 
   * method.
   *
   * The design of the Modifier class was inspired by the Lammps "Fix" 
   * base class, which provides a very flexible framework for designing
   * algorithms that can modify the state of a system.
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

      /**
      * Setup before entering the main loop.
      */ 
      virtual void setup(){};

      /// \name Integration action functions
      //@{ 
   
      virtual void preIntegrate1(long iStep) {};
      virtual void postIntegrate1(long iStep) {};
      virtual void preTransform(long iStep) {};
      virtual void preExchange(long iStep) {};
      virtual void postExchange(long iStep) {};
      virtual void postNeighbor(long iStep) {};
      virtual void preUpdate(long iStep) {};
      virtual void postUpdate(long iStep) {};
      virtual void preForce(long iStep) {};
      virtual void postForce(long iStep) {};
      virtual void endOfStep(long iStep) {};

      //@} 
      /// \name Interprocessor communication action functions
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

      /**
      * Return true if a flag is set, false otherwise.
      */
      bool isSet(Bit flag) const;

      /**
      * Return the unsigned int representation of all bit flags.
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

      /**
      * Call this to guarantee initialization of static variables.
      */
      static void initStatic();

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

      /// Unsigned int that stores a bitmap of all flags (one per bit).
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
