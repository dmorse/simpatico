#ifndef MCMD_ANALYZER_H
#define MCMD_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class
#include <util/misc/FileMaster.h>           // member variable

#include <string>
#include <iostream>
#include <fstream>

namespace McMd
{

   using namespace Util;

   /**
   * Abstract base for periodic output and/or analysis actions.
   *
   * The periodic action associated with an Analyzer can involve sampling
   * of a physical property and adding it to statistical accumulator, 
   * outputting it to file, or both. This periodic action must be 
   * implemented by the pure virtual sample() method.
   *
   * The sample() method should take the desired action only when the
   * simulation step index is an integer multiple of the associated interval
   * parameter.  The interval must be a positive integer that is a multiple 
   * of the static member Analyzer::baseInterval.
   *
   * The virtual sample() method does not take any parameters. An Analyzer
   * must thus access its parent Simulation and/or System via a pointer, 
   * which is usually initialized in its subclass constructor.
   *
   * Analyzer subclasses that are associated with one System, McSystem 
   * or MdSystem (i.e., almost all of them) should be derived from the
   * SystemAnalyzer<class SystemType> class template. This takes a 
   * reference to the parent system as parameter to its constructor. 
   * An Analyzer subclass that can be used with any System should be
   * derived from SystemAnalyzer<System>, and one that can be used 
   * only with a MdSystem or McSystem should be derived from 
   * SystemAnalyzer<MdSystem> or SystemAnalyzer<MdSystem>, 
   * respectively.
   *
   * \ingroup McMd_Analyzer_Module
   */
   class Analyzer : public ParamComposite
   {

   public:

      // Non-static Methods

      /**
      * Default constructor.
      */
      Analyzer();

      /**
      * Destructor.
      */
      virtual ~Analyzer();

      /**
      * Read parameters from archive.
      *
      * Default implementation, reads interval and outputFileName.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load parameters from archive.
      *
      * Default implementation, loads interval and outputFileName.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Load parameters to archive.
      *
      * Default implementation saves interval and outputFileName.
      *
      * \param ar loading/saving archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      *
      * Saves interval and outputFileName.
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Complete any required initialization.
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
      * The interval for an Analyzer must be a multiple of baseInterval.
      */
      static long baseInterval;

      /**
      * Define and initialize baseInterval.
      */
      static void initStatic();

   protected:

      /**
      * Set the FileMaster to use to open files.
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Read interval from file, with error checking.
      *
      * \param in input parameter file stream
      */
      void readInterval(std::istream &in);

      /**
      * Read outputFileName from file.
      *
      * \param in input parameter file stream.
      */
      void readOutputFileName(std::istream &in);

      /**
      * Load interval from archive, with error checking.
      *
      * \param ar input/loading archive
      */
      void loadInterval(Serializable::IArchive& ar);

      /**
      * Load output file name from archive.
      *
      * \param ar input/loading archive
      */
      void loadOutputFileName(Serializable::IArchive& ar);

      /**
      * Get the FileMaster by reference.
      *
      * This can be used to open multiple output files.
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

      /// Base name of output file(s).
      std::string outputFileName_;

      /// Number of simulation steps between subsequent actions.
      long interval_;

   private:

      /// Pointer to fileMaster for opening output file(s).
      FileMaster* fileMasterPtr_;

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

   // Method template

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void Analyzer::serialize(Archive& ar, const unsigned int version)
   {
      ar & interval_;
      ar & outputFileName_;
   }

}
#endif
