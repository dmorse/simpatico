#ifndef DDMD_SP_PROCESSOR_CPP
#define DDMD_SP_PROCESSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Processor.h"
#include <ddMd/sp/configIos/DdMdSpConfigIo.h>

namespace DdMd 
{

   /*
   * Constructor.
   */
   Processor::Processor()
    : configIoPtr_(0),
      configIoFactory_(*this),
      analyzerManager_(*this)
   {  setClassName("Processor"); }

   /*
   * Destructor.
   */
   Processor::~Processor()
   {
      if (configIoPtr_) {
         delete configIoPtr_;
      }
   }

   /*
   * Open, read and close a parameter file.
   */
   void Processor::readParam(const char* filename)
   {
      std::ifstream in;
      in.open(filename);
      readParam(in);
      in.close();
   }

   /*
   * Read parameters from file.
   */
   void Processor::readParameters(std::istream& in)
   {
      SpConfiguration::readParameters(in);
      readParamCompositeOptional(in, fileMaster_);
      readParamComposite(in, analyzerManager_);
   }

   bool Processor::hasFileMaster() const
   {  return fileMaster_.isActive(); }

   FileMaster& Processor::fileMaster()
   {
      if (!fileMaster_.isActive()) {
         UTIL_THROW(" Attempt to use inactive FileMaster");
      } 
      return fileMaster_;
   }

   // SpConfigIo management

   /*
   * Set SpConfigIo style.
   */
   void Processor::setSpConfigIo(std::string configIoName)
   {
      configIoPtr_ = configIoFactory_.factory(configIoName);
      if (configIoPtr_ == 0) {
         std::string msg;
         msg = "Unrecognized SpConfigIo subclass name: ";
         msg += configIoName;
         UTIL_THROW(msg.c_str());
      }
   }

   /**
   * Return the SpConfigIo (create default if necessary).
   */
   SpConfigIo& Processor::configIo() 
   {
      if (configIoPtr_ == 0) {
         configIoPtr_ = new DdMdSpConfigIo(*this);
         assert(configIoPtr_);
      }
      return *configIoPtr_;
   }

   // Reading and writing single configurations

   /*
   * Read a single configuration file.
   */
   void Processor::readConfig(std::ifstream& in)
   {
      clear();  
      configIo().readConfig(in); 
   }

   /*
   * Open, read and close a configuration file.
   */
   void Processor::readConfig(const char* filename)
   {
      std::ifstream inputFile;
      inputFile.open(filename);
      readConfig(inputFile);
      inputFile.close();
   }

   /*
   * Open, read and close a configuration file.
   */
   void Processor::readConfig(const std::string& filename)
   { readConfig(filename.c_str()); }

   /*
   * Write a single configuration file (must be open)
   */
   void Processor::writeConfig(std::ofstream& out)
   {  configIo().writeConfig(out); }

   /*
   * Open, write and close a configuration file.
   */
   void Processor::writeConfig(const std::string& filename)
   {
      std::ofstream outputFile;
      outputFile.open(filename.c_str());
      configIo().writeConfig(outputFile);
      outputFile.close();
   }

   // Analysis
   
   /*
   * Read and analyze a sequence of configuration files.
   */
   void Processor::analyzeDumps(int min, int max, std::string baseFileName)
   {
      // Preconditions
      if (min < 0)    UTIL_THROW("min < 0");
      if (max < min)  UTIL_THROW("max < min");

      //Timer timer;
      std::string filename;
      std::stringstream indexString;
      std::ifstream configFile;

      // Main loop
      Log::file() << "begin main loop" << std::endl;
      //timer.start();
      for (int iStep = min; iStep <= max; ++iStep) {

         indexString << iStep;
         filename = baseFileName;
         filename += indexString.str();
         configFile.open(filename.c_str());

         // Clear the stringstream
         indexString.str("");

         clear();
         readConfig(configFile);
         configFile.close();

         #if 0
         #ifndef INTER_NOPAIR
         // Build the configuration CellList
         pairPotential().buildCellList();
         #endif

         #ifdef UTIL_DEBUG
         isValid();
         #endif
         #endif

         // Initialize analyzers (taking in molecular information).
         if (iStep == min) {
            analyzerManager_.setup();
         }

         // Sample property values
         analyzerManager_.sample(iStep);

      }
      //timer.stop();
      Log::file() << "end main loop" << std::endl;

      // Output results of all analyzers to output files
      analyzerManager_.output();

      #if 0
      // Output time
      Log::file() << std::endl;
      Log::file() << "nConfig       " << nConfig << std::endl;
      Log::file() << "run time      " << timer.time()
                  << "  sec" << std::endl;
      Log::file() << "time / config " << timer.time()/double(nConfig)
                  << "  sec" << std::endl;
      Log::file() << std::endl;
      #endif

   }
}
#endif
