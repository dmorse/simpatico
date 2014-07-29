#ifndef DDMD_ANALYZER_H
#define DDMD_ANALYZER_H

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
   * The periodic action associated with an Analyzer is implemented 
   * by the pure virtual sample() method. This action often involves 
   * computation of a physical property which may be either added
   * added it to statistical accumulator or output to a file, or both. 
   *
   * The sample() method should take the desired action only when the
   * simulation step index is an integer multiple of the associated 
   * interval member variable. The interval must be a positive integer 
   * that is a multiple of the static member Analyzer::baseInterval.
   *
   * An Analyzer has access to its parent Simulation via the protected
   * Analyzer::simulation() method, which returns the parent Simulation 
   * by reference.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class Analyzer : public ParamComposite
   {

   public:

      // Non-static Methods

      /**
      * Constructor.
      */
      Analyzer(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~Analyzer();

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
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar)
      {}

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar)
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
      * The interval for an Analyzer must be a multiple of baseInterval.
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
      * is not a multiple of Analyzer::baseInterval, or if
      * baseInterval has not been set to a nonzero positive value.
      *
      * \param in input parameter file stream.
      */
      void readInterval(std::istream &in);

      /**
      * Load parameter interval from input archive.
      *
      * This function throws an exception if the value of interval
      * is not a multiple of Analyzer::baseInterval, or if
      * baseInterval has not been set to a nonzero positive value.
      *
      * \param ar input archive
      */
      void loadInterval(Serializable::IArchive &ar);

      /**
      * Save interval parameter to an archive.
      *
      * \param ar output archive
      */
      void saveInterval(Serializable::OArchive &ar);

      /**
      * Read outputFileName from file.
      *
      * \param in input parameter file
      */
      void readOutputFileName(std::istream &in);

      /**
      * Load output file name to an archive.
      *
      * \param ar input archive
      */
      void loadOutputFileName(Serializable::IArchive &ar);

      /**
      * Save output file name to an archive.
      *
      * \param ar output archive
      */
      void saveOutputFileName(Serializable::OArchive &ar);

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
   inline int Analyzer::interval() const
   {  return interval_; }

   /*
   * Return true iff the counter parameter is a multiple of the interval.
   */
   inline bool Analyzer::isAtInterval(long counter) const
   {  return (counter%interval_ == 0); }

   /*
   * Get the outputFileName string.
   */
   inline const std::string& Analyzer::outputFileName() const
   {  return outputFileName_; }

   /*
   * Get the parent Simulation by reference.
   */
   inline Simulation& Analyzer::simulation()
   {  return *simulationPtr_; }


}
#endif
