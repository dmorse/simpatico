#ifndef SPAN_PROCESSOR_CPP
#define SPAN_PROCESSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Processor.h"
#include <spAn/config/DdMdConfigReader.h>
#include <spAn/trajectory/LammpsDumpReader.h>

// std headers
#include <fstream>
#include <unistd.h>
#include <stdlib.h>

namespace SpAn 
{

   /*
   * Constructor.
   */
   Processor::Processor()
    : configReaderPtr_(0),
      trajectoryReaderPtr_(0),
      configReaderFactory_(*this),
      trajectoryReaderFactory_(*this),
      analyzerManager_(*this)
   {  setClassName("Processor"); }

   /*
   * Destructor.
   */
   Processor::~Processor()
   {
      if (configReaderPtr_) {
         delete configReaderPtr_;
      }
      if (trajectoryReaderPtr_) {
         delete trajectoryReaderPtr_;
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
      Configuration::readParameters(in);
      readParamCompositeOptional(in, fileMaster_);
      readParamComposite(in, analyzerManager_);
   }

   // ConfigReader Functions

   /*
   * Set ConfigReader style.
   */
   void Processor::setConfigReader(const std::string& configStyle)
   {
      configReaderPtr_ = configReaderFactory_.factory(configStyle);
      if (configReaderPtr_ == 0) {
         std::string msg;
         msg = "Unrecognized ConfigReader subclass name: ";
         msg += configStyle;
         UTIL_THROW(msg.c_str());
      }
   }

   // Config File Reader

   /*
   * Return the ConfigReader (create default if necessary).
   */
   ConfigReader& Processor::configReader() 
   {
      if (configReaderPtr_ == 0) {
         configReaderPtr_ = new DdMdConfigReader(*this);
         assert(configReaderPtr_);
      }
      return *configReaderPtr_;
   }

   /*
   * Read a single configuration file.
   */
   void Processor::readConfig(std::ifstream& in)
   {
      clear();  
      configReader().readConfig(in); 
   }

   /*
   * Open, read and close a configuration file.
   */
   void Processor::readConfig(const std::string& filename)
   { 
      std::ifstream inputFile;
      inputFile.open(filename.c_str());
      readConfig(inputFile);
      inputFile.close();
   }

   /*
   * Read and analyze a sequence of configuration files.
   */
   void Processor::analyzeDumps(int min, int max, const std::string& baseFileName)
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
         if (!configFile.is_open()) {
            std::string msg = "Configuration file is not open. Filename =";
            msg += filename;
            UTIL_THROW(msg.c_str());
         }

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

   // Trajectory Analysis

   /*
   * Set TrajectoryReader style.
   */
   void Processor::setTrajectoryReader(const std::string& trajectoryStyle)
   {
      trajectoryReaderPtr_ = 
                      trajectoryReaderFactory_.factory(trajectoryStyle);
      if (trajectoryReaderPtr_ == 0) {
         std::string msg;
         msg = "Unrecognized TrajectoryReader subclass name: ";
         msg += trajectoryStyle;
         UTIL_THROW(msg.c_str());
      }
   }

   /*
   * Return the TrajectoryReader (create default if necessary).
   */
   TrajectoryReader& Processor::trajectoryReader() 
   {
      if (trajectoryReaderPtr_ == 0) {
         trajectoryReaderPtr_ = new LammpsDumpReader(*this);
         assert(trajectoryReaderPtr_);
      }
      return *trajectoryReaderPtr_;
   }

   /*
   * Read and analyze a trajectory file.
   */
   void Processor::analyzeTrajectory(const std::string& filename)
   {

      // Open file
      std::ifstream file;
      file.open(filename.c_str());
      if (!file.is_open()) {
         std::string msg = "Trajectory file is not open. Filename =";
         msg += filename;
         UTIL_THROW(msg.c_str());
      }

      // Initialize analyzers (taking in molecular information).
      analyzerManager_.setup();

      // Main loop
      trajectoryReader().readHeader(file);
      Log::file() << "begin main loop" << std::endl;
      int iStep = 0;
      bool notEnd = true;
      while (notEnd) {
         notEnd = trajectoryReader().readFrame(file);
         if (notEnd) {
            analyzerManager_.sample(iStep);
            ++iStep;
         }
      }
      Log::file() << "end main loop" << std::endl;

      // Output any final results of analyzers to output files
      analyzerManager_.output();

      file.close();
   }

   /*
   * Does this processor have a FileMaster? 
   */
   bool Processor::hasFileMaster() const
   {  return fileMaster_.isActive(); }

   /*
   * Return FileMaster.
   */
   FileMaster& Processor::fileMaster()
   {
      if (!fileMaster_.isActive()) {
         UTIL_THROW(" Attempt to use inactive FileMaster");
      } 
      return fileMaster_;
   }

}
#endif
