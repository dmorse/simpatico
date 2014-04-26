#ifndef MDPP_PROCESSOR_CPP
#define MDPP_PROCESSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Processor.h"
#include "configIos/DdMdConfigIo.h"

namespace MdPp 
{

   /*
   * Constructor.
   */
   Processor::Processor()
    : configIoPtr_(0),
      configIoFactory_(*this),
      analyzerManager_(*this),
      newAtomPtr_(0)
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
   * Read parameters from file.
   */
   void Processor::readParameters(std::istream& in)
   {
      read<int>(in, "atomCapacity", atomCapacity_); 
      read<int>(in, "bondCapacity", bondCapacity_); 
      // etc. for angles, dihedrals

      atoms_.allocate(atomCapacity_);
      atomPtrs_.allocate(atomCapacity_);
      for (int i = 0; i < atomCapacity_; ++i) {
         atomPtrs_[i] = 0;
      }
      bonds_.allocate(bondCapacity_);
      // etc. for angles dihedrals

      readParamComposite(in, analyzerManager_);
   }

   // ConfigIo management

   /*
   * Set ConfigIo style.
   */
   void Processor::setConfigIo(std::string configIoName)
   {
      configIoPtr_ = configIoFactory_.factory(configIoName);
      if (configIoPtr_ == 0) {
         std::string msg;
         msg = "Unrecognized ConfigIo subclass name: ";
         msg += configIoName;
         UTIL_THROW(msg.c_str());
      }
   }

   /**
   * Return the ConfigIo (create default if necessary).
   */
   ConfigIo& Processor::configIo() 
   {
      if (configIoPtr_ == 0) {
         configIoPtr_ = new DdMdConfigIo(*this);
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
   void Processor::readConfig(const std::string& filename)
   {
      std::ifstream inputFile;
      inputFile.open(filename.c_str());
      readConfig(inputFile);
      inputFile.close();
   }

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

   // Mutators required by ConfigIos to read files.

   /*
   * Return pointer to location for new atom.
   */
   Atom* Processor::newAtomPtr()
   {
      if (newAtomPtr_) {
         UTIL_THROW("Error: an new atom is still active");
      }
      int size = atoms_.size() + 1;
      atoms_.resize(size);
      newAtomPtr_ = &atoms_[size - 1];
      return newAtomPtr_;
   }

   /*
   * Finalize addition of new atom.
   */
   void Processor::addAtom()
   {
      if (!newAtomPtr_) {
         UTIL_THROW("Error: No active new atom");
      }
      int id = newAtomPtr_->id;
      atomPtrs_[id] = newAtomPtr_;
      newAtomPtr_ = 0;
   }

   /*
   * Return pointer to location for new bond, and add to container.
   */
   Group<2>* Processor::newBondPtr()
   {
      int size = bonds_.size();
      bonds_.resize(size + 1);
      return &bonds_[size];
   }

   /*
   * Remove all atoms and bonds - set to empty state.
   */
   void Processor::clear()
   {
      atoms_.clear();
      bonds_.clear();
      for (int i = 0; i < atomCapacity_; ++i) {
         atomPtrs_[i] = 0;
      }
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
         // Build the system CellList
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
