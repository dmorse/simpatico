#ifndef MDPP_ANALYZER_H
#define MDPP_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/FileMaster.h>       // member variable

#include <string>
#include <iostream>
#include <fstream>

namespace Util { class FileMaster; }

namespace MdPp
{

   using namespace Util;
   class Configuration;
   class Processor;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * The periodic action of an Analyzer can involve sampling of a
   * physical property and adding it to statistical accumulator, 
   * outputting it to file, or both. This periodic action must be 
   * implemented by the pure virtual sample() method.
   *
   * The sample() method should take the desired action only when the
   * configuration step index is an integer multiple of the associated 
   * interval parameter.  
   *
   * The virtual sample() method does not take any parameters. An 
   * Analyzer must thus access its parent Configuration via a pointer, 
   * which is usually initialized in its subclass constructor.
   *
   * \ingroup MdPp_Analyzer_Module
   */
   class Analyzer : public ParamComposite
   {

   public:

      // Non-static Methods

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      */
      Analyzer(Configuration& configuration);

      /**
      * Constructor.
      *
      * \param processor parent Processor object
      */
      Analyzer(Processor& processor);

      /**
      * Constructor.
      *
      * \param configuration parent Configuration
      * \param fileMaster associated Util::FileMaster object
      */
      Analyzer(Configuration& configuration, FileMaster& fileMaster);

      /**
      * Destructor.
      */
      virtual ~Analyzer();

      /**
      * Setup before configuration.
      *
      * This method is called just before the beginning of the
      * main loop. It may be used to complete any initialization 
      * or checks that require knowledge of the configuration. 
      *
      * The default implementation is empty.
      */
      virtual void setup()
      {}

      /**
      * Clear statistical accumulators.
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
      * \param iStep current step index.
      */
      virtual void sample(long iStep) = 0;

      /**
      * Output any results at the end of an analysis.
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
      * \param counter configuration step counter
      */
      bool isAtInterval(long counter) const;

      // Static members

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
      * Read outputFileName from file.
      *
      * \param in input parameter file
      */
      void readOutputFileName(std::istream &in);

      /**
      * Get the parent Configuration by reference.
      */
      Configuration& configuration();

      /**
      * Get an associated FileMaster by reference.
      */
      FileMaster& fileMaster();

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

      /// Pointer to parent Configuration
      Configuration* configurationPtr_;

      /// Pointer to FileMaster
      FileMaster* fileMasterPtr_;

      /// Number of configuration steps between subsequent actions.
      long interval_;

      /// True if this is reponsible for FileMaster destruction.
      bool ownsFileMaster_;

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
   * Get the parent Configuration by reference.
   */
   inline Configuration& Analyzer::configuration()
   {  return *configurationPtr_; }

}
#endif
