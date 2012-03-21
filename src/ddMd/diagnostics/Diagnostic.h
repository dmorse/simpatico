#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class
#include <ddMd/util/FileMaster.h>           // member variable

#include <string>
#include <iostream>
#include <fstream>

namespace DdMd
{

   using namespace Util;
   class System;

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
   * must thus access its parent Simulation and/or System via a pointer, 
   * which is usually initialized in its subclass constructor.
   *
   * Diagnostic subclasses that are associated with one System, McSystem 
   * or MdSystem (i.e., almost all of them) should be derived from the
   * SystemDiagnostic<class SystemType> class template. This takes a 
   * reference to the parent system as parameter to its constructor. 
   * A Diagnostic subclass that can be used with any System should be
   * derived from SystemDiagnostic<System>, and one that can be used 
   * only with a MdSystem or McSystem should be derived from 
   * SystemDiagnostic<MdSystem> or SystemDiagnostic<MdSystem>, 
   * respectively.
   *
   * \ingroup Diagnostic_Module
   */
   class Diagnostic : public ParamComposite
   {

   public:

      // Non-static Methods

      /**
      * Default constructor.
      */
      Diagnostic(System& system);

      /**
      * Destructor.
      */
      virtual ~Diagnostic();

      /**
      * Setup before simulation.
      *
      * This method must be called just before the beginning of
      * the main simulation loop, after an initial configuration 
      * is known. It may be used to complete any initialization
      * that cannot be completed in the readParam method, because
      * knowledge of the configuration is needed. 
      *
      * The default implementation is an empty function.
      */
      virtual void setup()
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
      * Output any results at the end of the simulation.
      *
      * The default implementation is an empty function.
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
      * Get the parent System by reference.
      */
      System& system();

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

      /// Pointer to parent System
      System* systemPtr_;

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
   * Get the parent System by reference.
   */
   inline System& Diagnostic::system()
   {  return *systemPtr_; }



}
#endif
