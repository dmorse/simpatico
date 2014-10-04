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
#include <spAn/config/DdMdConfigWriter.h>
#include <spAn/trajectory/LammpsDumpReader.h>
#include <util/format/Str.h>

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
      configWriterFactory_(*this),
      trajectoryReaderFactory_(*this),
      analyzerManager_(*this),
      paramFileName_(),
      commandFileName_()
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
   * Process command line options.
   */
   void Processor::setOptions(int argc, char * const * argv)
   {
      bool  eFlag = false;               // option e

      // Read command-line arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "e")) != -1) {
         switch (c) {
         case 'e':
           eFlag = true;
           break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
         }
      }

      // If option -e, enable echoing of parameters as they are read
      if (eFlag) {
         ParamComponent::setEcho(true);
      }

      // Check number of remaining arguments
      int nArg = argc - optind;
      if (nArg < 1) {
         Log::file() << "optind = " << optind << std::endl;
         Log::file() << "argc   = " << argc << std::endl;
         UTIL_THROW("Two few non-option arguments");
      }
      if (nArg > 2) {
         Log::file() << "optind = " << optind << std::endl;
         Log::file() << "argc   = " << argc << std::endl;
         UTIL_THROW("Two many non-option arguments");
      }
   
      // Read parameter file name
      paramFileName_ = argv[optind];
      Log::file() << "paramFileName  = " << paramFileName_ << std::endl;

      // Read command file name, if any
      if (nArg == 2) {
         commandFileName_ = argv[optind+1];
         Log::file() << "commandFileName = " 
                     << commandFileName_ << std::endl;
      }

   }

   /*
   * Read default param file.
   */
   void Processor::readParam()
   {
      if (paramFileName_ != "") {
         readParam(paramFileName_.c_str());
      } else {
         UTIL_THROW("No param file specified");
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

   /*
   * Read and execute commands from default command file.
   */
   void Processor::readCommands()
   {
      if (commandFileName_ != "") {
         std::ifstream file;
         fileMaster().openInputFile(commandFileName_, file);
         readCommands(file);
         file.close();
      } else {
         readCommands(std::cin);
      }
   }

   /*
   * Read and execute commands from a specified command file.
   */
   void Processor::readCommands(std::istream &in)
   {
      #if 0
      if (!isInitialized_) {
         UTIL_THROW("Processor is not initialized");
      }
      #endif

      std::string command;
      std::string filename;
      std::ifstream inputFile;
      std::ofstream outputFile;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "SET_CONFIG_READER") {
            std::string classname;
            in >> classname;
            Log::file() << " " << classname << std::endl;
            setConfigReader(classname);
         } else
         if (command == "READ_CONFIG") {

            // If needed, read auxiliary file
            if (configReader().needsAuxiliaryFile()) {
               Log::file() << "\nReading auxiliary file: "; 
               in >> filename;
               Log::file() << filename;
               fileMaster().openInputFile(filename, inputFile);
               configReader().readAuxiliaryFile(inputFile);
               inputFile.close();
            }

            // Read actual configuration file
            Log::file() << "\nReading config file: ";
            in >> filename;
            Log::file() << filename << std::endl;
            fileMaster().openInputFile(filename, inputFile);
            configReader().readConfig(inputFile);
            inputFile.close();

         } else
         if (command == "SET_CONFIG_WRITER") {
            std::string classname;
            in >> classname;
            Log::file() << " " << classname << std::endl;
            setConfigWriter(classname);
         } else
         if (command == "WRITE_CONFIG") {
            if (configWriter().needsAuxiliaryFile()) {
               Log::file() << "\nReading auxiliary file: ";
               in >> filename;
               Log::file() << filename;
               fileMaster().openInputFile(filename, inputFile);
               configWriter().readAuxiliaryFile(inputFile);
               inputFile.close();
            }
            Log::file() << "\nWriting config file: ";
            in >> filename;
            Log::file() << filename << std::endl;
            fileMaster().openOutputFile(filename, outputFile);
            configWriter().writeConfig(outputFile);
            outputFile.close();
         } else
         if (command == "SET_TRAJECTORY_READER") {
            std::string classname;
            in >> classname;
            Log::file() << " " << classname << std::endl;
            setTrajectoryReader(classname);
         } else
         if (command == "ANALYZE_TRAJECTORY") {
            in >> filename;
            Log::file() << " " << filename << std::endl;
            analyzeTrajectory(filename);
         } else 
         {
            Log::file() << "  Error: Unknown command  " << std::endl;
            readNext = false;
         }

      }
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
      // Log::file() << "Setting config reader to class "
      //            << configReaderPtr_->className() << std::endl;
   }

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
      std::ifstream file;
      file.open(filename.c_str());
      readConfig(file);
      file.close();
   }

   /*
   * Read and analyze a sequence of configuration files.
   */
   void Processor::analyzeConfigs(const std::string& baseFileName,
                                  int min, int max, int interval)
   {
      // Preconditions
      if (min < 0)    UTIL_THROW("min < 0");
      if (max < min)  UTIL_THROW("max < min");
      if (interval <= 0)  UTIL_THROW("interval <= 0");

      //Timer timer;
      std::string filename;
      std::stringstream indexString;
      std::ifstream configFile;

      // Main loop
      Log::file() << "begin main loop" << std::endl;
      //timer.start();
      for (int iStep = min; iStep <= max; iStep += interval) {

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

   // ConfigWriter Functions

   /*
   * Set ConfigWriter style.
   */
   void Processor::setConfigWriter(const std::string& configStyle)
   {
      configWriterPtr_ = configWriterFactory_.factory(configStyle);
      if (configWriterPtr_ == 0) {
         std::string msg;
         msg = "Unrecognized ConfigWriter subclass name: ";
         msg += configStyle;
         UTIL_THROW(msg.c_str());
      }
   }

   /*
   * Return the ConfigWriter (create default if necessary).
   */
   ConfigWriter& Processor::configWriter() 
   {
      if (configWriterPtr_ == 0) {
         configWriterPtr_ = new DdMdConfigWriter(*this);
         assert(configWriterPtr_);
      }
      return *configWriterPtr_;
   }

   /*
   * Read a configuration file.
   */
   void Processor::writeConfig(std::ofstream& in)
   {  configWriter().writeConfig(in); }

   /*
   * Open, write and close a configuration file.
   */
   void Processor::writeConfig(const std::string& filename)
   { 
      std::ofstream file;
      file.open(filename.c_str());
      writeConfig(file);
      file.close();
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
      if (trajectoryReader().isBinary()) {
         fileMaster().openInputFile(filename, file,
                                    std::ios::in | std::ios::binary);
      } else {
         fileMaster().openInputFile(filename, file);
      }
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
      //if (!fileMaster_.isActive()) {
      //   UTIL_THROW(" Attempt to use inactive FileMaster");
      //} 
      return fileMaster_;
   }

}
#endif
