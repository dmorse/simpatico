#ifndef MDPP_PROCESSOR_H
#define MDPP_PROCESSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/storage/Configuration.h>                // base class
#include <mdPp/config/ConfigReaderFactory.h>           // member 
#include <mdPp/config/ConfigWriterFactory.h>           // member 
#include <mdPp/trajectory/TrajectoryReaderFactory.h>   // member 
#include <mdPp/processor/ProcessorAnalyzerManager.h>   // member 
#include <util/misc/FileMaster.h>                       // member 

namespace MdPp 
{

   class ConfigReader;

   using namespace Util;

   /**
   * A post-processor for analyzing outputs of MD simulations.
   *
   * \ingroup MdPp_Storage_Module
   */
   class Processor : public Configuration
   {

   public:

      /**
      * Constructor
      */
      Processor();

      /**
      * Destructor
      */
      ~Processor();

      /// \name Initialization
      //@{
      
      /**
      * Process command line options.
      *  
      * \param argc number of arguments
      * \param argv array of argument C-strings
      */
      void setOptions(int argc, char * const * argv);

      using ParamComposite::readParam;

      /**
      * Read parameter file specified in command line, if any.
      *
      * Does nothing if no parameter file has been specified.
      */
      void readParam();

      /**
      * Open, read, and close parameter file.
      *
      * \param filename input parameter file name
      */
      void readParam(const char* filename);

      /**
      * Read parameters.
      *
      * \param in input command file stream
      */
      void readParameters(std::istream& in);

      /**
      * Read command file specified in command line.
      * 
      * If no command file has been specified, reads from the command line.
      */
      void readCommands();

      /**
      * Read and execute commands from a command file.
      *
      * \param in input command file stream
      */
      void readCommands(std::istream &in);

      //@}
      /// \name ConfigReader Interface 
      //@{
      
      /**
      * Set ConfigReader style  (creates a ConfigReader).
      *
      * \param configReaderName identifier for ConfigReader subclass
      */
      void setConfigReader(const std::string& configReaderName);

      /**
      * Return the current ConfigReader (create default if necessary).
      */
      ConfigReader& configReader();
   
      /**
      * Read a single configuration file.
      */
      void readConfig(std::ifstream& in);

      /**
      * Open, read and close a configuration file.
      */
      void readConfig(const std::string& filename);
   
      /**
      * Read and analyze a sequence of numbered configuration files.
      *
      * This function reads and analyzes a sequence of configuration files 
      * that were generated by running a previous simulation. The function 
      * reads files with names of the form inputPrefix() + n for integer 
      * suffixes min <= n <= max with subsequent values differing by the
      * specified interval.
      *
      * \param baseFileName  root name for dump files (without int suffix)
      * \param min  integer suffix of first configuration file name
      * \param max  integer suffix of last configuration file name
      * \param interval  interval between subsequent timestep values
      */  
      void analyzeConfigs(const std::string& baseFileName,
                          int min, int max, int interval = 1);

      //@}
      /// \name ConfigWriter Interface 
      //@{
      
      /**
      * Set ConfigWriter style  (creates a ConfigWriter).
      *
      * \param configWriterName identifier for ConfigWriter subclass
      */
      void setConfigWriter(const std::string& configWriterName);

      /**
      * Return the current ConfigWriter (create default if necessary).
      */
      ConfigWriter& configWriter();

      /**
      * Write a single configuration file.
      */
      void writeConfig(std::ofstream& in);

      /**
      * Open, write and close a configuration file.
      */
      void writeConfig(const std::string& filename);
   
      //@}
      /// \name Trajectory File Interface
      //@{

      /**
      * Set TrajectoryReader style  (creates a TrajectoryReader).
      *
      * \param trajectoryStyle TrajectoryReader subclass identifier
      */
      void setTrajectoryReader(const std::string& trajectoryStyle);

      /**
      * Return the current TrajectoryReader (create default if necessary).
      */
      TrajectoryReader& trajectoryReader();
   
      /**
      * Open, read, analyze and close a single trajectory file.
      *
      * \param filename name of trajectory file.
      */
      void analyzeTrajectory(const std::string& filename);

      //@}
      /// \name Miscellaneous functions
      //@{
 
      /**
      * Return FileMaster if active, or throw Exception.
      */
      FileMaster& fileMaster();

      //@}

   private:

      /// Pointer to current ConfigReader object.
      ConfigReader* configReaderPtr_;

      /// Pointer to current ConfigWriter object.
      ConfigWriter* configWriterPtr_;

      /// Pointer to current TrajectoryReader object.
      TrajectoryReader* trajectoryReaderPtr_;

      /// Factory for choosing ConfigReader at run time.
      ConfigReaderFactory configReaderFactory_;

      /// Factory for choosing ConfigWriter at run time.
      ConfigWriterFactory configWriterFactory_;

      /// Factory for choosing TrajectoryReader at run time.
      TrajectoryReaderFactory trajectoryReaderFactory_;

      /// Manager for analyzers
      ProcessorAnalyzerManager analyzerManager_;

      /// FileMaster
      FileMaster fileMaster_;

      /// String identifier for ConfigReader class name
      std::string configReaderName_;

   };

}
#endif
