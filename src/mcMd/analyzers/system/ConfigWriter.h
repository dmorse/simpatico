#ifndef MCMD_CONFIG_WRITER_H
#define MCMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <mcMd/simulation/System.h>

namespace McMd
{

   using namespace Util;

   /**
   * Periodically dump system configuration to a new file.
   *
   * This class uses the System::writeConfig() method to write the entire
   * system configuration periodically, at interval read from file. Names 
   * of configuration dump files are constructed by concatenating the
   * System outputPrefix, a dumpPrefix, the literal string "config.", and 
   * an integer that is incremented after each dump. The dump prefix is a 
   * string that is read by the read() method of this class. For example, 
   * if System::outputPrefix is "out/" and dumpPrefix is "dump/", the 
   * configurations will be written to files "out/dump/config.0",
   * "out/dump/config.1", etc.
   *
   * \sa \ref mcMd_analyzer_ConfigWriter_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_Module
   */
   class ConfigWriter : public SystemAnalyzer<System>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System object. 
      */
      ConfigWriter(System& system);
   
      /**
      * Destructor.
      */
      virtual ~ConfigWriter()
      {} 
   
      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear nSample counter.
      */
      virtual void setup();
  
      /**
      * Dump configuration to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);
  
   private:
      
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long nSample_;
   
      /// Has readParam been called?
      long isInitialized_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void ConfigWriter::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & nSample_;
   }

}
#endif 
