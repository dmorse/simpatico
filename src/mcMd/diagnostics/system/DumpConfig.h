#ifndef MCMD_DUMP_CONFIG_H
#define MCMD_DUMP_CONFIG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>
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
   * \ingroup McMd_Diagnostic_Module
   */
   class DumpConfig : public SystemDiagnostic<System>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System object. 
      */
      DumpConfig(System& system);
   
      /**
      * Destructor.
      */
      virtual ~DumpConfig()
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
      long    nSample_;
   
      /// Has readParam been called?
      long    isInitialized_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void DumpConfig::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSample_;
   }

}
#endif 
