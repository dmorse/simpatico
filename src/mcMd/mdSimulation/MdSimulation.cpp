#ifndef MD_SIMULATION_CPP
#define MD_SIMULATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdDiagnosticManager.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/trajectoryIos/TrajectoryIo.h>
#include <mcMd/species/Species.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/util/Log.h>
#include <util/random/serialize.h>
#include <util/archives/Serializable_includes.h>
#include <util/util/ioUtil.h>
#include <util/util/Timer.h>


#include <sstream>
#include <iostream>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MdSimulation::MdSimulation()
    : Simulation(),
      system_(),
      mdDiagnosticManagerPtr_(0),
      isRestarting_(false)
   {
      system_.setId(0);
      system_.setSimulation(*this);
      system_.setFileMaster(fileMaster());
      mdDiagnosticManagerPtr_ = new MdDiagnosticManager(*this);
      assert(mdDiagnosticManagerPtr_);
      setDiagnosticManager(mdDiagnosticManagerPtr_);
   }

   /* 
   * Destructor.
   */
   MdSimulation::~MdSimulation()
   {
      assert(mdDiagnosticManagerPtr_);
      if (mdDiagnosticManagerPtr_) {
         delete mdDiagnosticManagerPtr_;
         mdDiagnosticManagerPtr_ = 0;
      }
   }

   /* 
   * Read parameters from file.
   */
   void MdSimulation::readParam(std::istream &in)
   { 
      readBegin(in, "MdSimulation");

      Simulation::readParam(in); 

      readParamComposite(in, system_); 
      readParamComposite(in, diagnosticManager());

      isValid();

      readEnd(in);
   }
 
   /*
   * Read default parameter file.
   */
   void MdSimulation::readParam()
   {  readParam(fileMaster().paramFile()); }


   /*
   * Read and implement commands in an input script.
   */
   void MdSimulation::readCommands(std::istream &in)
   {
      std::string   command;
      std::string   filename;
      std::ifstream inputFile;
      std::ofstream outputFile;

      #ifndef UTIL_MPI
      std::istream&     inBuffer = in;
      #else
      std::stringstream inBuffer;
      std::string       line;
      #endif

      bool readNext = true;
      while (readNext) {

         #ifdef UTIL_MPI
         if (!hasParamCommunicator() || isParamIoProcessor()) {
            getNextLine(in, line);
         } 
         if (hasParamCommunicator()) {
            bcast<std::string>(communicator(), line, 0);
         }
         inBuffer.clear();
         for (unsigned i=0; i < line.size(); ++i) {
            inBuffer.put(line[i]);
         }
         #endif

         inBuffer >> command;
         Log::file().setf(std::ios_base::left);
         Log::file().width(15);
         Log::file() << command;

         if (isRestarting_) {

            if (command == "RESTART") {
               int endStep;
               inBuffer >> endStep;
               Log::file() << "  " << iStep_ << " to " 
                           << endStep << std::endl;
               simulate(endStep, isRestarting_);
               isRestarting_ = false;
            } else {
               UTIL_THROW("Missing RESTART command");
            }

         } else {

            if (command == "READ_CONFIG") {
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openInputFile(filename, inputFile);
               system().readConfig(inputFile);
               inputFile.close();
            } else
            if (command == "THERMALIZE") {
               double temperature;
               inBuffer >> temperature;
               Log::file() << Dbl(temperature, 15, 6) << std::endl;
               system().setBoltzmannVelocities(temperature);
               system().removeDriftVelocity();
            } else
            if (command == "SIMULATE") {
               int endStep;
               inBuffer >> endStep;
               Log::file() << Int(endStep, 15) << std::endl;
               bool isContinuation = false;
               simulate(endStep, isContinuation);
            } else
            if (command == "CONTINUE") {
               if (iStep_ == 0) {
                  UTIL_THROW("Attempt to continue simulation when iStep_ == 0");
               }
               int endStep;
               inBuffer >> endStep;
               Log::file() << Int(endStep, 15) << std::endl;
               bool isContinuation = true;
               simulate(endStep, isContinuation);
            } else
            if (command == "ANALYZE_DUMPS") {
               int min, max;
               inBuffer >> min >> max >> filename;
               Log::file() <<  Int(min, 15) << Int(max, 15)
                           <<  Str(filename, 20) << std::endl;
               analyzeDumps(min, max, filename);
            } else
            if (command == "WRITE_CONFIG") {
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openOutputFile(filename, outputFile);
               system().writeConfig(outputFile);
               outputFile.close();
            } else
            if (command == "WRITE_PARAM") {
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openOutputFile(filename, outputFile);
               writeParam(outputFile);
               outputFile.close();
            } else
            if (command == "SET_CONFIG_IO") {
               std::string classname;
               inBuffer >> classname;
               Log::file() << Str(classname, 15) << std::endl;
               system().setConfigIo(classname);
            } else
            if (command == "ANALYZE_TRAJECTORY") {
               std::string classname;
               std::string filename;
               int min, max;
               inBuffer >> min >> max >> classname >> filename;
               Log::file() << " " << Str(classname,15) 
                           << " " << Str(filename, 15)
                           << std::endl;
               analyzeTrajectory(min, max, classname,filename);
            } else 
            if (command == "GENERATE_MOLECULES") {
               double boxL;
               DArray<double> exclusionRadius;
               DArray<int> nMolecule;
               exclusionRadius.allocate(nAtomType());
               nMolecule.allocate(nSpecies());
               inBuffer >> boxL;
               Log::file() << "  " << boxL;
               for (int iSpecies=0; iSpecies < nSpecies(); iSpecies++) {
                  inBuffer >> nMolecule[iSpecies];
                  Log::file() << "  " << nMolecule[iSpecies];
               }
               for (int iType=0; iType < nAtomType(); iType++) {
                  inBuffer >> exclusionRadius[iType];
                  Log::file() << "  " << exclusionRadius[iType];
               }
               Log::file() << std::endl;
               Boundary boundary;
               boundary.setCubicLengths(boxL);
               for (int iSpecies=0; iSpecies < nSpecies(); iSpecies++) {
                  species(iSpecies).generateMolecules(
                     nMolecule[iSpecies], exclusionRadius, system(),
                     &system().bondPotential(),
                     boundary);
               }
            } else
            if (command == "FINISH") {
               Log::file() << std::endl;
               readNext = false;
            } else {
               Log::file() << "Error: Unknown command  " << std::endl;
               readNext = false;
            }

         }
      }
   }

   /*
   * Read and implement commands from the default command file.
   */
   void MdSimulation::readCommands()
   {  readCommands(fileMaster().commandFile()); }

   /* 
   * Run a simulation.
   */
   void MdSimulation::simulate(int endStep, bool isContinuation)
   {
      Timer timer;

      if (isContinuation) {
         Log::file() << "Continuing simulation from iStep = " 
                     << iStep_ << std::endl;
      } else {
         iStep_ = 0;
         diagnosticManager().initialize();
         system_.mdIntegrator().setup();
      }
      int beginStep = iStep_;
      int nStep = endStep - beginStep;

      // Main loop 
      Log::file() << std::endl;
      timer.start();
      for ( ; iStep_ < endStep; ++iStep_) {

         // Sample scheduled diagnostics
         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               system().shiftAtoms();
               diagnosticManager().sample(iStep_);
            }
         }

         #ifdef MCMD_NOPAIR
         else {
            // When the pair potential is disabled, require that
            // Diagnostic::baseInterval != 0 to guarantee periodic 
            // shifting of atomic positions into primary cell.
            UTIL_THROW("Diagnostic::baseInterval == 0 with no pair potential");
         }
         #endif

         // Take one MD step with the MdIntegrator
         system_.mdIntegrator().step();

      }
      timer.stop();
      double time  = timer.time();
      double rstep = double(nStep);

      // Shift final atomic positions 
      system().shiftAtoms();

      // Sample diagnostic for final configuration
      if (Diagnostic::baseInterval > 0) {
         if (iStep_ % Diagnostic::baseInterval == 0) {
            diagnosticManager().sample(iStep_);
         }
      }

      // Final diagnostic output
      diagnosticManager().output();

      // Output timing results
      Log::file() << std::endl;
      Log::file() << "Time Statistics" << std::endl;
      Log::file() << "endStep              " << endStep << std::endl;
      Log::file() << "nStep                " << nStep   << std::endl;
      Log::file() << "run time             " 
                  << time << " sec" << std::endl;
      Log::file() << "time / nStep         " 
                  << time / rstep << " sec" << std::endl;
      Log::file() << "time / (nStep*nAtom) " 
                  <<  time / (rstep*double(system().nAtom())) 
                  << " sec" << std::endl;
      Log::file() << std::endl;
      Log::file() << std::endl;

      #ifndef MCMD_NOPAIR
      Log::file() << "PairList Statistics" << std::endl;
      Log::file() << "maxNPair           " 
                  << system().pairPotential().pairList().maxNPair()
                  << std::endl;
      Log::file() << "maxNAtom           " 
                  << system().pairPotential().pairList().maxNAtom()
                  << std::endl;
      int buildCounter = system().pairPotential().pairList().buildCounter();
      Log::file() << "buildCounter       " 
                  << buildCounter
                  << std::endl;
      Log::file() << "steps / build      " 
                  << rstep/double(buildCounter)
                  << std::endl;
      Log::file() << std::endl;
      #endif

   }

   /*
   * Read and analyze a sequence of configuration files.
   */
   void MdSimulation::analyzeDumps(int min, int max, std::string dumpPrefix)
   {
      // Preconditions
      if (min < 0)    UTIL_THROW("min < 0");
      if (max < min)  UTIL_THROW("max < min");

      Timer             timer;
      std::string       filename;
      std::stringstream indexString;
      std::ifstream     configFile;
      int               nConfig;
      nConfig = max - min + 1;

      // Main loop
      Log::file() << "begin main loop" << std::endl;
      timer.start();
      for (iStep_ = min; iStep_ <= max; ++iStep_) {

         indexString << iStep_;
         filename = dumpPrefix;
         filename += indexString.str();
         fileMaster().openInputFile(filename, configFile);

         // Clear the stringstream
         indexString.str("");
        
         system().readConfig(configFile);
         configFile.close();

         #ifndef MCMD_NOPAIR
         // Build the system CellList
         system().pairPotential().buildPairList();
         #endif

         #ifndef NDEBUG
         isValid();
         #endif

         // Initialize diagnostics (taking in molecular information).
         if (iStep_ == min) diagnosticManager().initialize();

         // Sample property values
         diagnosticManager().sample(iStep_);

         // Clear out the System for the next readConfig.
         system().removeAllMolecules();
 
      }
      timer.stop();
      Log::file()<< "end main loop" << std::endl;

      // Output results of all diagnostics to output files
      diagnosticManager().output();

      // Output time 
      Log::file() << std::endl;
      Log::file() << "nConfig       " << nConfig << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / config " << timer.time()/double(nConfig) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;

   }

   /*
   * Read and analyze a trajectory file
   */
   void MdSimulation::analyzeTrajectory(int min, int max, std::string classname,
      std::string filename)
   {
      Timer             timer;
      std::stringstream indexString;
      std::fstream*     trajectoryFile;
      TrajectoryIo*     trajectoryIo;
      int nFrames; 

      // Obtain trajectoryIo objecct
      if (!(trajectoryIo = system().trajectoryIoFactory().factory(classname))) {
        std::string message;
        message="Invalid TrajectoryIo class name " + classname;
        UTIL_THROW(message.c_str());
      }

      Log::file() << "reading " << filename << std::endl;

      // Open trajectory file
      trajectoryFile = new std::fstream();
      trajectoryFile->open(filename.c_str(), std::ios::in | std::ios::binary);

      if (trajectoryFile->fail()) {
        std::string message;
        message= "Error opening trajectory file. Filename: " + filename;
        UTIL_THROW(message.c_str());
      }

      // read in information from trajectory file
      trajectoryIo->readHeader(*trajectoryFile);

      nFrames = trajectoryIo->nFrames();

      // Main loop
      Log::file() << "begin main loop" << std::endl;
      timer.start();

      // Allow for negative values of min, max (counted from the end)
      if (min < 0) min = nFrames+min;
      if (max < 0) max = nFrames+max;

      // Sanity checks
      if (min < 0 || min >= nFrames)  UTIL_THROW("min < 0 or min >= nFrames!");
      if (max < 0 || max >= nFrames)  UTIL_THROW("max < 0 or max >= nFrames!");
      if (max < min)  UTIL_THROW("max < min!");

      for (iStep_ = 0; iStep_ <= max; ++iStep_) {

         // Read frames, even if they are not sampled
         trajectoryIo->readFrame(*trajectoryFile);

         #ifndef MCMD_NOPAIR
         // Build the system CellList
         system().pairPotential().buildPairList();
         #endif

         #ifndef NDEBUG
         isValid();
         #endif

         // Initialize diagnostics (taking in molecular information).
         if (iStep_ == min) diagnosticManager().initialize();

         // Sample property values only for iStep >= min
         if (iStep_ >= min) diagnosticManager().sample(iStep_);

      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;

      // Close trajectory file
      trajectoryFile->close();

      // delete objects
      delete trajectoryIo;
      delete trajectoryFile;

      // Output results of all diagnostics to output files
      diagnosticManager().output();

      // Output time 
      Log::file() << std::endl;
      Log::file() << "nFrames       " << max-min+1 << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;

   }

   void MdSimulation::writeRestart(const std::string& filename)
   {
      std::ofstream out;
      fileMaster().openParamOFile(filename, ".prm", out);
      writeParam(out);
      out << "iStep = " << iStep_; 
      out.close();

      fileMaster().openRestartOFile(filename, ".rst", out);
      Serializable::OArchiveType ar;
      ar.setStream(out);

      ar & random();
      //ar & system();
      system().save(ar);
      system().mdIntegrator().save(ar);
      ar & iStep_;
      mdDiagnosticManagerPtr_->save(ar);
      out.close();
   }

   void MdSimulation::readRestart(const std::string& filename)
   {
      isRestarting_ = true;

      // Open and read parameter (*.prm) file
      std::ifstream in;
      fileMaster().openParamIFile(filename, ".prm", in);
      readParam(in);
      in.close();

      // Open restart *.rst file and associate with an archive
      fileMaster().openRestartIFile(filename, ".rst", in);
      Serializable::IArchiveType ar;
      ar.setStream(in);

      // Load internal state from restart (*.rst) file.
      ar & random();
      system().load(ar);
      system().mdIntegrator().load(ar);
      ar & iStep_;
      mdDiagnosticManagerPtr_->load(ar);
      in.close();

      // Read command (*.cmd) file
      fileMaster().openParamIFile(filename, ".cmd", in);
      readCommands(in);
      in.close();
   }

   /* 
   * Return true if this MdSimulation is valid, or throw an Exception.
   */
   bool MdSimulation::isValid() const
   {
      Simulation::isValid();
      system_.isValid();

      if (&system_.simulation() != this) {
         UTIL_THROW("Invalid pointer to Simulation in System");
      }

      return true;
   }

}    
#endif
