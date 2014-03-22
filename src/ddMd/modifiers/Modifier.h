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
   * A Modifier can modify the time evolution of the simulation.
   *
   * Modifier is an abstract base class for classes that can modify the
   * system during an MD simulation, and thus modify the time evolution
   * defined by the integrator. The class declares a variety of virtual
   * methods that, when re-implemented and activated by subclasses, will
   * be invoked at specific points during the main integration loop, thus
   * allowing the designer of a subclass to add almost arbitrary actions
   * to the underlying integration algorithm. 
   *
   * The Modifier class declares interfaces for several virtual functions
   * that that the integrator can be instructed to invoke at specific
   * points within the body of the main integration loop. We refer to 
   * theses as integrator action functions. The name of each such function 
   * describes when it will be invoked. For example, the postForce() method 
   * is invoked (if at all) immediately after evaluation of all forces 
   * within the main integration loop.  Modifier also declares several 
   * virtual functions with names that contain the words "pack" and "unpack",
   * which we will call communication functions, that may be defined by
   * subclasses to pack additional information into the messages that 
   * are communicated between processors during exchange and update steps.
   *
   * Each subclass of Modifier will normally implement only a few of 
   * these possible action and communication methods. To record which
   * virtual functions it defines, each subclass must also "activate" 
   * each such function by invoking the set() method within its 
   * constructor passing this function a named constant that activates
   * a particular virtual method.  Thus, for example, a subclass that 
   * implements the virtual postForce() function must also call the
   * function set(Flags::PostForce) in its constructor. The constant 
   * Flags::PostForce and other named constants are defined as static 
   * constant members of the Modifier::Flags nested class. The name
   * of each such constant is given by a capitalized version of the
   * name of the corresponding virtual method. When the set() function
   * is invoked to activate an integrator or communication function,
   * it adds that function to list of functions that will be invoked 
   * at the appropriate point in the integration algorithm. The 
   * constants associated with communcation functions each activate 
   * a pair of associated "pack" and "unpack" methods that must work
   * together pack additional information into an MPI message on one
   * processor and unpack it at the receiving processor. 
   *
   * Each Modifier also has an integer member named interval that
   * determines the interval, in time steps, with which integrator 
   * actions should be taken. An activated integrator action function 
   * will be executed only during time steps that are multiples of 
   * this interval. Each Modifier object has a single interval value,
   * which applies to all action methods that it implements and
   * activates. The value of the interval is initialized to 1 (i.e.,
   * every time step) in the base Modifier class constructor, but may 
   * be reset in the subclass readParameters method by invoking the 
   * readInterval() method to read a value from file. 
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
 
      /** 
      * Call just before the first step of velocity-Verlet algorithm. 
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void preIntegrate1(long iStep) {};

      /** 
      * Call just after the first step of velocity-Verlet algorithm. 
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void postIntegrate1(long iStep) {};

      /** 
      * Call on exchange steps before transforming to scaled atomic coordinates.
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void preTransform(long iStep) {};

      /** 
      * Call on exchange steps after transforming but before exchanging atoms.
      *
      * Atom positions are scaled [0,1] on entry and return.
      */
      virtual void preExchange(long iStep) {};

      /** 
      * Call on exchange steps right after atom exchange, but before reneighboring
      *
      * Atom positions are scaled [0,1] on entry and return.
      */
      virtual void postExchange(long iStep) {};

      /** 
      * Call on exchange steps after re-building neighbor list (reneighboring).
      *
      * Atom positions are scaled [0,1] on entry and return.
      */
      virtual void postNeighbor(long iStep) {};

      /** 
      * Call on update steps before updating ghost positions.
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void preUpdate(long iStep) {};

      /** 
      * Call on update steps after updating ghost positions.
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void postUpdate(long iStep) {};

      /** 
      * Call after updating but before calculating forces.
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void preForce(long iStep) {};

      /** 
      * Call after calculating forces
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void postForce(long iStep) {};

      /** 
      * Call at the end of the time step.
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void endOfStep(long iStep) {};

      //@} 
      /// \name Interprocessor communication action functions
      //@{ 

      /**
      * Pack data into buffer used to exchange atoms.
      */
      virtual void packExchange() {};

      /**
      * Unpack data from buffer used to exchange atoms.
      */
      virtual void unpackExchange() {};

      /**
      * Pack data into buffer used to update ghost positions.
      */
      virtual void packUpdate() {};

      /**
      * Unpack data from buffer used to update ghost positions.
      */
      virtual void unpackUpdate() {};

      /**
      * Pack data into buffer used to reverse update forces.
      *
      * Will only be used if reverse communication is enabled.
      */
      virtual void packReverseUpdate() {};

      /**
      * Unpack data from the buffer used to reverse update forces.
      *
      * Will only be used if reverse communication is enabled.
      */
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

         /// Flag to activate setup() function.
         static const Bit Setup;

         /// Flag to activate preIntegrate() function.
         static const Bit PreIntegrate1;

         /// Flag to activate postIntegrate1() function.
         static const Bit PostIntegrate1;

         /// Flag to activate preTransform() function.
         static const Bit PreTransform;

         /// Flag to activate preExchange() function.
         static const Bit PreExchange;

         /// Flag to activate postExchange() function.
         static const Bit PostExchange;

         /// Flag to activate postNeighbor() function.
         static const Bit PostNeighbor;

         /// Flag to activate preUpdate() function.
         static const Bit PreUpdate;

         /// Flag to activate postUpdate() function.
         static const Bit PostUpdate;

         /// Flag to activate preForce() function.
         static const Bit PreForce;

         /// Flag to activate postForce() function.
         static const Bit PostForce;

         /// Flag to activate endOfStep() function.
         static const Bit EndOfStep;

         /// Flag to activate pack/unpack exchange functions.
         static const Bit Exchange;

         /// Flag to activate pack/unpack update functions.
         static const Bit Update;

         /// Flag to activate pack/unpack reverse update functions.
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
