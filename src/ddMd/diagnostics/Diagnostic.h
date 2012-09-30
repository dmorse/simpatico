#ifndef DDMD_DIAGNOSTIC_H
#define DDMD_DIAGNOSTIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class
#include <util/misc/FileMaster.h>           // member variable

#include <string>
#include <iostream>
#include <fstream>

namespace DdMd
{

   using namespace Util;
   class Simulation;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * The periodic action associated with a Diagnostic can involve sampling
   * of a physical property and adding it to statistical accumulator, 
   * outputting it to file, or both. This periodic action must be 
   * implemented by the pure virtual sample() method.
   *
   * The sample() method should take the desired action only when the
   * simulation step index is an integer multiple of the associated interval
   * parameter.  The interval must be a positive integer that is a multiple 
   * of the static member Diagnostic::baseInterval.
   *
   * The virtual sample() method does not take any parameters. A Diagnostic
   * must thus access its parent Simulation and/or Simulation via a pointer, 
   * which is usually initialized in its subclass constructor.
   *
   * \ingroup DdMd_Diagnostic_Module
   */
   class Diagnostic : public ParamComposite
   {

   public:

      // Non-static Methods

      /**
      * Default constructor.
      */
      Diagnostic(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Diagnostic();

      /**
      * Setup before simulation.
      *
      * This method is called just before the beginning of the
      * main simulation loop within the Integrator::run() method.
      * It may be used to complete any initialization or checks 
      * that require knowledge of the configuration. It will be 
      * called every time run is invoked, not just the first.
      *
      * The default implementation is empty.
      */
      virtual void setup()
      {}

      /**
      * Clear statistical accumulators.
      *
      * This method is called by the Integrator::clear() method,
      * which is called before the main simulation loop only the 
      * first time that the integrator is run, and when invoked
      * explicitly thereafter.
      *
      * The default implementation is empty.
      */
      virtual void clear()
      {}

      /**
      * Calculate, analyze and/or output a physical quantity.
      *
      * Take an action if iStep is a multiple of interval.
      * If iStep is not a multiple of interval, this method
      * should do nothing and return immediately.
      *
      * \param iStep current simulation step index.
      */
      virtual void sample(long iStep) = 0;

      /**
      * Output any results at the end of a simulation.
      *
      * The default implementation is empty.
      */
      virtual void output()
      {}

      /**
      * Get interval value.
      */
      int interval() const;

      /**
      * Return true iff counter is a multiple of the interval.
      *
      * \param counter simulation step counter
      */
      bool isAtInterval(long counter) const;

      // Static members

      /**
      * The interval for a Diagnostic must be a multiple of baseInterval.
      */
      static long baseInterval;

      /**
      * Define and initialize baseInterval.
      */
      static void initStatic();

   protected:

      /**
      * Read parameter interval from file.
      *
      * This function throws an exception if the value of interval
      * is not a multiple of Diagnostic::baseInterval, or if
      * baseInterval has not been set to a nonzero positive value.
      *
      * \param in input parameter file stream.
      */
      void readInterval(std::istream &in);

      /**
      * Read outputFileName from file.
      *
      * \param in input parameter file stream.
      */
      void readOutputFileName(std::istream &in);

      /**
      * Get the parent Simulation by reference.
      */
      Simulation& simulation();

      /**
      * Return outputFileName string.
      */
      const std::string& outputFileName() const;

      /**
      * Return outputFileName string with added suffix.
      */
      std::string outputFileName(const std::string& suffix) const;

   private:

      /// Base name of output file(s).
      std::string outputFileName_;

      /// Pointer to parent Simulation
      Simulation* simulationPtr_;

      /// Number of simulation steps between subsequent actions.
      long   interval_;

   };

   // Inline methods

   /*
   * Return interval value.
   */
   inline int Diagnostic::interval() const
   {  return interval_; }

   /*
   * Return true iff the counter parameter is a multiple of the interval.
   */
   inline bool Diagnostic::isAtInterval(long counter) const
   {  return (counter%interval_ == 0); }

   /*
   * Get the outputFileName string.
   */
   inline const std::string& Diagnostic::outputFileName() const
   {  return outputFileName_; }

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& Diagnostic::simulation()
   {  return *simulationPtr_; }


}
#endif
