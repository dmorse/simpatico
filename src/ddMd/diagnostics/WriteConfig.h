#ifndef DDMD_WRITE_CONFIG_H
#define DDMD_WRITE_CONFIG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/Diagnostic.h>
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically dump simulation configuration to a new file.
   *
   * This class uses the Simulation::writeConfig() method to write the entire
   * simulation configuration periodically, at interval read from file. Names 
   * of configuration dump files are constructed by concatenating the
   * Simulation outputPrefix, a dumpPrefix, the literal string "config.", and 
   * an integer that is incremented after each dump. The dump prefix is a 
   * string that is read by the read() method of this class. For example, 
   * if Simulation::outputPrefix is "out/" and dumpPrefix is "dump/", the 
   * configurations will be written to files "out/dump/config.0",
   * "out/dump/config.1", etc.
   *
   * \ingroup DdMd_Diagnostic_Module
   */
   class WriteConfig : public Diagnostic
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      WriteConfig(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~WriteConfig()
      {} 
   
      /**
      * Read parameters and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
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

}
#endif 
