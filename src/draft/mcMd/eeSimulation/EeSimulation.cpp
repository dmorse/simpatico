/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "EeSimulation.h"
#include <mcMd/mcSimulation/McAnalyzerManager.h>
#include <mcMd/simulation/serialize.h>
#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/mcMoves/McMoveManager.h>
#include <mcMd/species/Species.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#include <mcMd/perturb/ReplicaMove.h>
#endif
#endif

#include <util/param/Factory.h>
#include <util/random/serialize.h>
#include <util/archives/Serializable_includes.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/format/Str.h>
#include <util/misc/Log.h>
#include <util/misc/Timer.h>
#include <util/misc/ioUtil.h>

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace McMd
{

   using namespace Util;

   #ifdef UTIL_MPI
   /*
   * Constructor.
   */
   EeSimulation::EeSimulation(MPI::Intracomm& communicator)
    : Simulation(communicator),
      system_(),
      mcMoveManagerPtr_(0),
      mcAnalyzerManagerPtr_(0),
      paramFilePtr_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      // Set connections between this EeSimulation and child McSystem
      system().setId(0);
      system().setSimulation(*this);
      system().setFileMaster(fileMaster());

      // Create McMove and Analyzer managers
      mcMoveManagerPtr_ = new McMoveManager(*this);
      mcAnalyzerManagerPtr_ = new McAnalyzerManager(*this);

      // Pass Manager<Analyzer>* to Simulation base class.
      setAnalyzerManager(mcAnalyzerManagerPtr_);
   }
   #endif

   /*
   * Constructor.
   */
   EeSimulation::EeSimulation()
    : Simulation(),
      system_(),
      mcMoveManagerPtr_(0),
      mcAnalyzerManagerPtr_(0),
      paramFilePtr_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      // Set connections between this EeSimulation and child McSystem
      system().setId(0);
      system().setSimulation(*this);
      system().setFileMaster(fileMaster());

      // Create McMove and Analyzer managers
      mcMoveManagerPtr_ = new McMoveManager(*this);
      mcAnalyzerManagerPtr_ = new McAnalyzerManager(*this);

      // Pass Manager<Analyzer>* to Simulation base class.
      setAnalyzerManager(mcAnalyzerManagerPtr_);
   }

   /*
   * Destructor.
   */
   EeSimulation::~EeSimulation()
   {
      delete mcMoveManagerPtr_;
      delete mcAnalyzerManagerPtr_;
   }

   /*
   * Read parameters from file.
   */
   void EeSimulation::readParam(std::istream &in)
   {
      // Record identity of parameter file
      paramFilePtr_ = &in;

      readBegin(in,"EeSimulation");

      // Read all species, analyzers, random number seed
      Simulation::readParam(in);

      // Read the McSystem parameters: potential parameters, temperature etc.
      readParamComposite(in, system());

      // Read Monte Carlo Moves
      assert(mcMoveManagerPtr_);
      readParamComposite(in, *mcMoveManagerPtr_);

      // Read Analyzers
      readParamComposite(in, analyzerManager());

      isValid();
      isInitialized_ = true;
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   void EeSimulation::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read and execute commands from a specified command file.
   */
   void EeSimulation::readCommands(std::istream &in)
   {
      if (!isInitialized_) {
         UTIL_THROW("EeSimulation is not initialized");
      }

      std::string    command;
      std::string    filename;
      std::ifstream  inputFile;
      std::ofstream  outputFile;

      #ifndef UTIL_MPI
      std::istream&     inBuffer = in;
      #else
      std::stringstream inBuffer;
      std::string       line;
      #endif

      bool readNext = true;
      while (readNext) {

         // fdef   UTIL_MPI, read a line and copy to stringstream inBuffer
         // ifndef UTIL_MPI, inBuffer is simply a reference to istream in.

         #ifdef UTIL_MPI
         // Read a command line, and broadcast if necessary.
         if (!hasIoCommunicator() || isIoProcessor()) {
            getNextLine(in, line);
         }
         if (hasIoCommunicator()) {
            bcast<std::string>(communicator(), line, 0);
         }

         // Copy the command line to inBuffer.
         inBuffer.clear();
         for (unsigned i=0; i < line.size(); ++i) {
            inBuffer.put(line[i]);
         }
         #endif

         inBuffer >> command;
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

            if (command == "FINISH") {
               Log::file() << std::endl;
               readNext = false;
            } else
            if (command == "SET_CONFIG_IO") {
               std::string classname;
               inBuffer >> classname;
               Log::file() << Str(classname, 15) << std::endl;
               system().setConfigIo(classname);
            } else
            if (command == "READ_CONFIG") {
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openInputFile(filename, inputFile);
               system().readConfig(inputFile);
               inputFile.close();
            } else
            if (command == "SIMULATE") {
               int endStep;
               inBuffer >> endStep;
               Log::file() << "  " << endStep << std::endl;
               bool isContinuation = false;
               simulate(endStep, isContinuation);
            } else
            if (command == "CONTINUE") {
               if (iStep_ == 0) {
                  UTIL_THROW("Attempt to continue when iStep_ == 0");
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
               Log::file() << "  " <<  min << "  " <<  max
                           << "  " <<  filename << std::endl;
               analyze(min, max, filename);
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
               Log::file() << "  " << filename << std::endl;
               fileMaster().openOutputFile(filename, outputFile);
               writeParam(outputFile);
               outputFile.close();
            } else 
            if (command == "GENERATE_MOLECULES") {
               double boxL;
               DArray<double> exclusionRadius;
               DArray<int> nMolecule;
               exclusionRadius.allocate(nAtomType());
               nMolecule.allocate(nSpecies());
               inBuffer >> boxL;
               Log::file() << "  " << boxL;
               for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
                  inBuffer >> nMolecule[iSpecies];
                  Log::file() << "  " << nMolecule[iSpecies];
               }
               for (int iType=0; iType < nAtomType(); iType++) {
                  inBuffer >> exclusionRadius[iType];
                  Log::file() << "  " << exclusionRadius[iType];
               }
               Log::file() << std::endl;
               Boundary boundary;
               boundary.setCubic(boxL);
               for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
                  species(iSpecies).generateMolecules(
                     nMolecule[iSpecies], exclusionRadius, system(),
                     &system().bondPotential(),
                     boundary);   
               }
               #ifndef SIMP_NOPAIR 
               // Generate cell list
               system().pairPotential().buildCellList();
               #endif
            } else {
               Log::file() << "  Error: Unknown command  " << std::endl;
               readNext = false;
            }

         }
      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   void EeSimulation::readCommands()
   {  readCommands(fileMaster().commandFile()); }

   /*
   * Run this MC simulation.
   */
   void EeSimulation::simulate(int endStep, bool isContinuation)
   {
      if (!isInitialized_) {
         UTIL_THROW("EeSimulation not initialized");
      }

      if (isContinuation) {
         Log::file() << "Restarting from iStep = " << iStep_ << std::endl;
      } else {
         iStep_ = 0;
         analyzerManager().setup();
         mcMoveManagerPtr_->setup();
      }
      int beginStep = iStep_;
      int nStep = endStep - beginStep;
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for ( ; iStep_ < endStep; ++iStep_) {

         // Sample analyzers
         if (Analyzer::baseInterval > 0) {
            if (iStep_ % Analyzer::baseInterval == 0) {
               analyzerManager().sample(iStep_);
            }
         }

         // Choose and attempt an McMove
         mcMoveManagerPtr_->chooseMove().move();

         #ifdef MCMD_PERTURB
         #ifdef UTIL_MPI
         // Attempt replica move, if any.
         if (system().hasPerturbation()) {
            if (system().hasReplicaMove()) {
               if (system().replicaMove().isAtInterval(iStep_)) {
                  bool success = system().replicaMove().move();
                  #ifndef SIMP_NOPAIR
                  if (success) {
                     system().pairPotential().buildCellList();
                  }
                  #endif
               }
            }
         }
         #endif
         #endif

      }
      timer.stop();
      double time = timer.time();

      // Final analyzer sample
      assert(iStep_ == endStep);
      if (Analyzer::baseInterval > 0) {
         if (iStep_ % Analyzer::baseInterval == 0) {
            analyzerManager().sample(iStep_);
         }
      }

      // Output results of all analyzers to output files
      analyzerManager().output();

      // Output results of move statistics to files
      mcMoveManagerPtr_->output();

      // Output time for the run
      Log::file() << std::endl;
      Log::file() << "endStep       " << endStep << std::endl;
      Log::file() << "nStep         " << nStep << std::endl;
      Log::file() << "run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep  " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Print McMove acceptance statistics
      long attempt;
      long accept;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(32) << left <<  "Move Name"
           << setw(12) << right << "Attempted"
           << setw(12) << right << "Accepted"
           << setw(15) << right << "AcceptRate"
           << endl;
      int nMove = mcMoveManagerPtr_->size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = (*mcMoveManagerPtr_)[iMove].nAttempt();
         accept  = (*mcMoveManagerPtr_)[iMove].nAccept();
         Log::file() << setw(32) << left 
              << mcMoveManagerPtr_->className(iMove)
              << setw(12) << right << attempt
              << setw(12) << accept
              << setw(15) << fixed << setprecision(6)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

      #ifdef MCMD_PERTURB
      #ifdef UTIL_MPI
      // Print replica-exchange acceptance statistics
      if (system().hasPerturbation()) {
         if (system().hasReplicaMove()) {
  
            double ratio; 
            int    nAttempt, nAccept;
   
            nAttempt = system().replicaMove().nAttempt();
            nAccept = system().replicaMove().nAccept();
            ratio = nAttempt == 0 ? 0.0 : double(nAccept)/double(nAttempt);
            Log::file() << "Replica Exchange "
                        << Int(nAttempt) << Int(nAccept) << Dbl(ratio)
                        << std::endl;
         }
      }
      #endif
      #endif

   }

   /*
   * Read and analyze a sequence of configuration files.
   */
   void EeSimulation::analyze(int min, int max, std::string dumpPrefix)
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

         #ifndef SIMP_NOPAIR
         // Build the system CellList
         system().pairPotential().buildCellList();
         #endif

         #ifdef UTIL_DEBUG
         isValid();
         #endif

         // Initialize analyzers (taking in molecular information).
         if (iStep_ == min) analyzerManager().setup();

         // Sample property values
         analyzerManager().sample(iStep_);

         // Clear out the System for the next readConfig.
         system().removeAllMolecules();

      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;

      // Output results of all analyzers to output files
      analyzerManager().output();

      // Output time
      Log::file() << std::endl;
      Log::file() << "nConfig       " << nConfig << std::endl;
      Log::file() << "run time      " << timer.time()
                  << "  sec" << std::endl;
      Log::file() << "time / config " << timer.time()/double(nConfig)
                  << "  sec" << std::endl;
      Log::file() << std::endl;

      iStep_ = 0;
   }

   void EeSimulation::writeRestart(const std::string& filename)
   {
      std::ofstream out;
      fileMaster().openParamOFile(filename, ".prm", out);
      writeParam(out);
      out.close();

      fileMaster().openRestartOFile(filename, ".rst", out);
      Serializable::OArchive ar;
      ar.setStream(out);
      ar & random();
      ar & system();
      ar & iStep_;
      mcMoveManagerPtr_->save(ar);
      mcAnalyzerManagerPtr_->save(ar);
      out.close();
   }

   void EeSimulation::readRestart(const std::string& filename)
   {
      isRestarting_ = true;
      std::ifstream in;

      // Open and read parameter (*.prm) file
      fileMaster().openParamIFile(filename, ".prm", in);
      readParam(in);
      in.close();

      // Open restart (*.rst) file and associate with an archive
      fileMaster().openRestartIFile(filename, ".rst", in);
      Serializable::IArchive ar;
      ar.setStream(in);

      // Load state from restart file
      ar & random();
      system().load(ar);
      ar & iStep_;
      mcMoveManagerPtr_->load(ar);
      mcAnalyzerManagerPtr_->load(ar);
      in.close();

      fileMaster().openParamIFile(filename, ".cmd", in);
      readCommands(in);
      in.close();

   }

   /*
   * Get the McMove factory.
   */
   Factory<McMove>& EeSimulation::mcMoveFactory()
   {  return mcMoveManagerPtr_->factory(); }

   /*
   * Check validity: return true if valid, or throw Exception.
   */
   bool EeSimulation::isValid() const
   {
      Simulation::isValid();
      system().isValid();

      if (&system().simulation() != this) {
         UTIL_THROW("Invalid pointer to Simulation in System");
      }

      return true;
   }

}
