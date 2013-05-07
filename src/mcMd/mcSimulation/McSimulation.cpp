#ifndef MCMD_MC_SIMULATION_CPP
#define MCMD_MC_SIMULATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "McSimulation.h"
#include "McDiagnosticManager.h"
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/mcMoves/McMoveManager.h>
#include <mcMd/species/Species.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef UTIL_MPI
#ifdef MCMD_PERTURB
#include <mcMd/perturb/ReplicaMove.h>
#endif
#endif

#include <util/param/Factory.h>
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
#include <unistd.h>

namespace McMd
{

   using namespace Util;

   #ifdef UTIL_MPI
   /*
   * Constructor.
   */
   McSimulation::McSimulation(MPI::Intracomm& communicator)
    : Simulation(communicator),
      system_(),
      mcMoveManagerPtr_(0),
      mcDiagnosticManagerPtr_(0),
      paramFilePtr_(0),
      saveFileName_(),
      saveInterval_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("McSimulation"); 

      // Set connections between this McSimulation and child McSystem
      system().setId(0);
      system().setSimulation(*this);
      system().setFileMaster(fileMaster());

      // Create McMove and Diagnostic managers
      mcMoveManagerPtr_ = new McMoveManager(*this);
      mcDiagnosticManagerPtr_ = new McDiagnosticManager(*this);

      // Pass Manager<Diagnostic>* to Simulation base class.
      setDiagnosticManager(mcDiagnosticManagerPtr_);
   }
   #endif

   /*
   * Constructor.
   */
   McSimulation::McSimulation()
    : Simulation(),
      system_(),
      mcMoveManagerPtr_(0),
      mcDiagnosticManagerPtr_(0),
      paramFilePtr_(0),
      saveFileName_(),
      saveInterval_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("McSimulation"); 

      // Set connections between this McSimulation and child McSystem
      system().setId(0);
      system().setSimulation(*this);
      system().setFileMaster(fileMaster());

      // Create McMove and Diagnostic managers
      mcMoveManagerPtr_ = new McMoveManager(*this);
      mcDiagnosticManagerPtr_ = new McDiagnosticManager(*this);

      // Pass Manager<Diagnostic>* to Simulation base class.
      setDiagnosticManager(mcDiagnosticManagerPtr_);
   }

   /*
   * Destructor.
   */
   McSimulation::~McSimulation()
   {
      delete mcMoveManagerPtr_;
      delete mcDiagnosticManagerPtr_;
   }

   /*
   * Process command line options.
   */
   void McSimulation::setOptions(int argc, char **argv)
   {
      char* rarg   = 0;
      bool  eflag  = false;
      bool  rflag  = false;
      #ifdef MCMD_PERTURB
      bool  pflag = false;
      #endif
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "epr:")) != -1) {
         switch (c) {
         case 'e':
           eflag = true;
           break;
         case 'r':
           rflag = true;
           rarg  = optarg;
           break;
         #ifdef MCMD_PERTURB
         case 'p':
           pflag = true;
           break;
         #endif
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }
   
      // Set flag to echo parameters as they are read.
      if (eflag) {
         Util::ParamComponent::setEcho(true);
      }

      #ifdef MCMD_PERTURB
      // Set to use a perturbation.
      if (pflag) {

         if (rflag) {
            std::string msg("Error: Options -r and option -p are incompatible. Use -r alone. ");
            msg += "Existence of a perturbation is specified in restart file.";
            UTIL_THROW(msg.c_str());
         }
   
         // Set to expect perturbation in the param file.
         system().setExpectPerturbation();
   
         #ifdef UTIL_MPI
         Util::Log::file() << "Set to use single parameter and command files" 
                           << std::endl;
         setIoCommunicator();
         #endif
   
      }
      #endif

      // Read a restart file.
      if (rflag) {
         Log::file() << "Begin reading restart, base file name " 
                     << std::string(rarg) << std::endl;
         isRestarting_ = true; 
         load(std::string(rarg));
         Util::Log::file() << std::endl;
      }

   }

   /*
   * Read parameters from file.
   */
   void McSimulation::readParameters(std::istream &in)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called readParam when already initialized");
      }

      // Record identity of parameter file
      paramFilePtr_ = &in;

      // Read all species, diagnostics, random number seed
      Simulation::readParameters(in);

      // Read the McSystem parameters: potential parameters, temperature etc.
      readParamComposite(in, system());

      // Read Monte Carlo Moves
      assert(mcMoveManagerPtr_);
      readParamComposite(in, *mcMoveManagerPtr_);

      // Read Diagnostics
      readParamComposite(in, diagnosticManager());

      // Parameters for writing restart files
      read<int>(in, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         read<std::string>(in, "saveFileName", saveFileName_);
      }

      isValid();
      isInitialized_ = true;
   }

   /*
   * Read parameter block, including begin and end.
   */
   void McSimulation::readParam(std::istream& in)
   {
      if (isRestarting_) {
         if (isInitialized_) {
            return;
         }
      }
      readBegin(in, className().c_str());
      readParameters(in);
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   void McSimulation::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Load internal state from an archive.
   */
   void McSimulation::loadParameters(Serializable::IArchive &ar)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called readParam when already initialized");
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         UTIL_THROW("Error: Has a param communicator in loadParameters");
      }
      #endif

      Simulation::loadParameters(ar);
      loadParamComposite(ar, system());
      loadParamComposite(ar, *mcMoveManagerPtr_);
      loadParamComposite(ar, diagnosticManager());
      loadParameter<int>(ar, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         loadParameter<std::string>(ar, "saveFileName", saveFileName_);
      }

      system().loadConfig(ar);
      ar >> iStep_;
      isValid();
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void McSimulation::save(Serializable::OArchive &ar)
   {
      Simulation::save(ar);
      system().saveParameters(ar);
      mcMoveManagerPtr_->save(ar);
      diagnosticManager().save(ar);
      ar << saveInterval_;
      if (saveInterval_ > 0) {
         ar << saveFileName_;
      }

      system().saveConfig(ar);
      ar << iStep_;
   }

   void McSimulation::load(const std::string& filename)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called load when already initialized");
      }
      if (!isRestarting_) {
         UTIL_THROW("Error: Called load without restart option");
      }

      // Load from archive
      Serializable::IArchive ar;
      fileMaster().openRestartIFile(filename, ".rst", ar.file());
      load(ar);
      ar.file().close();

      // Set command (*.cmd) file
      std::string commandFileName = filename + ".cmd";
      fileMaster().setCommandFileName(commandFileName);

      #ifdef UTIL_MPI
      #ifdef MCMD_PERTURB
      if (system().hasPerturbation()) {
         // Read one command file, after reading multiple restart files.
         Util::Log::file() << "Set to use a single command file" 
                           << std::endl;
         setIoCommunicator();
      }
      #endif
      #endif

      isInitialized_ = true;
   }

   void McSimulation::save(const std::string& filename)
   {
      if (saveInterval_ > 0) {
         if (iStep_ % saveInterval_ == 0) {
            Serializable::OArchive ar;
            fileMaster().openRestartOFile(filename, ".rst", ar.file());
            save(ar);
            ar.file().close();
         }
      }
   }

   /*
   * Read and execute commands from a specified command file.
   */
   void McSimulation::readCommands(std::istream &in)
   {
      if (!isInitialized_) {
         UTIL_THROW("McSimulation is not initialized");
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

         // iffdef UTIL_MPI, read a line and copy to stringstream inBuffer
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
               DArray<double> ExclusionRadii;
               DArray<int>    capacities;
               ExclusionRadii.allocate(nAtomType());
               capacities.allocate(nSpecies());

               // Parse command
               inBuffer >> system().boundary();
               Log::file() << "  " << system().boundary();
               Label capacityLabel("Capacities:");
               inBuffer >> capacityLabel;
               for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
                  inBuffer >> capacities[iSpecies];
                  Log::file() << "  " << capacities[iSpecies];
               }
               Label radiusLabel("Radii:");
               inBuffer >> radiusLabel;
               for (int iType=0; iType < nAtomType(); iType++) {
                  inBuffer >> ExclusionRadii[iType];
                  Log::file() << "  " << ExclusionRadii[iType];
               }
               Log::file() << std::endl;

               for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
                  species(iSpecies).generateMolecules(
                     capacities[iSpecies], ExclusionRadii, system(),
                     &system().bondPotential(),
                     system().boundary());   
               }

               #ifndef INTER_NOPAIR 
               // Generate cell list
               system().pairPotential().buildCellList();
               #endif

            } else
            if (command == "DEFORM_CELL") {
               
               // Read in configuration from file
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openInputFile(filename, inputFile);
               system().readConfig(inputFile);
               inputFile.close();

               System::MoleculeIterator molIter;
               Molecule::AtomIterator atomIter;
               for (int iSpec=0; iSpec < nSpecies(); ++iSpec) {
                  system().begin(iSpec, molIter);
                  for ( ; molIter.notEnd(); ++molIter) {
                     molIter->begin(atomIter);
                     for ( ; atomIter.notEnd(); ++atomIter) {
                        Vector cartPosition, genPosition;
                        cartPosition = atomIter->position();
                        system().boundary().transformCartToGen(cartPosition, genPosition);
                        atomIter->position() = genPosition;
                     }
                  }
               }

               // Read in new boundary
               inBuffer >> system().boundary();
               Log::file() << "  " << system().boundary();
               Log::file() << std::endl;

               for (int iSpec=0; iSpec < nSpecies(); ++iSpec) {
                  system().begin(iSpec, molIter);
                  for ( ; molIter.notEnd(); ++molIter) {
                     molIter->begin(atomIter);
                     for ( ; atomIter.notEnd(); ++atomIter) {
                        Vector cartPosition, genPosition;
                        genPosition = atomIter->position();
                        system().boundary().transformGenToCart(genPosition, cartPosition);
                        atomIter->position() = cartPosition;
                     }
                  }
               }

               // Write out configuration to file
               inBuffer >> filename;
               Log::file() << Str(filename, 15) << std::endl;
               fileMaster().openOutputFile(filename, outputFile);
               system().writeConfig(outputFile);
               outputFile.close();

               #ifndef INTER_NOPAIR 
               // Generate cell list
               system().pairPotential().buildCellList();
               #endif

            } else
            #ifndef UTIL_MPI
            #ifndef INTER_NOPAIR
            if (command == "SET_PAIR") {
               std::string paramName;
               int typeId1, typeId2; 
               double value;
               inBuffer >> paramName >> typeId1 >> typeId2 >> value;
               Log::file() << "  " <<  paramName 
                           << "  " <<  typeId1 << "  " <<  typeId2
                           << "  " <<  value << std::endl;
               system().pairPotential()
                       .set(paramName, typeId1, typeId2, value);
            } else 
            #endif // ifndef INTER_NOPAIR
            if (command == "SET_BOND") {
               std::string paramName;
               int typeId; 
               double value;
               inBuffer >> paramName >> typeId >> value;
               Log::file() << "  " <<  paramName << "  " <<  typeId 
                           << "  " <<  value << std::endl;
               system().bondPotential().set(paramName, typeId, value);
            } else 
            #ifdef INTER_ANGLE
            if (command == "SET_ANGLE") {
               std::string paramName;
               int typeId; 
               double value;
               inBuffer >> paramName >> typeId >> value;
               Log::file() << "  " <<  paramName << "  " <<  typeId 
                           << "  " <<  value << std::endl;
               system().anglePotential().set(paramName, typeId, value);
            } else 
            #endif // ifdef INTER_ANGLE
            #ifdef INTER_DIHEDRAL
            if (command == "SET_DIHEDRAL") {
               std::string paramName;
               int typeId; 
               double value;
               inBuffer >> paramName >> typeId >> value;
               Log::file() << "  " <<  paramName << "  " <<  typeId 
                           << "  " <<  value << std::endl;
               system().dihedralPotential().set(paramName, typeId, value);
            } else 
            #endif // ifdef INTER_DIHEDRAL
            #endif // ifndef UTIL_MPI
            {
               Log::file() << "  Error: Unknown command  " << std::endl;
               readNext = false;
            }

         }
      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   void McSimulation::readCommands()
   {  readCommands(fileMaster().commandFile()); }

   /*
   * Run this MC simulation.
   */
   void McSimulation::simulate(int endStep, bool isContinuation)
   {
      if (!isInitialized_) {
         UTIL_THROW("McSimulation not initialized");
      }

      if (isContinuation) {
         Log::file() << "Restarting from iStep = " 
                     << iStep_ << std::endl;
      } else {
         iStep_ = 0;
         diagnosticManager().setup();
         mcMoveManagerPtr_->setup();
      }
      int beginStep = iStep_;
      int nStep = endStep - beginStep;
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for ( ; iStep_ < endStep; ++iStep_) {

         // Diagnostics and restart outut
         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               save(saveFileName_);
               diagnosticManager().sample(iStep_);
            }
         }

         // Choose and attempt an McMove
         mcMoveManagerPtr_->chooseMove().move();

         #ifdef UTIL_MPI
         #ifdef MCMD_PERTURB
         // Attempt replica move, if any.
         if (system().hasPerturbation()) {
            if (system().hasReplicaMove()) {
               if (system().replicaMove().isAtInterval(iStep_)) {
                  bool success = system().replicaMove().move();
                  #ifndef INTER_NOPAIR
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

      // Final diagnostics / save
      assert(iStep_ == endStep);
      if (Diagnostic::baseInterval > 0) {
         if (iStep_ % Diagnostic::baseInterval == 0) {
            save(saveFileName_);
            diagnosticManager().sample(iStep_);
         }
      }

      // Output results of all diagnostics to output files
      diagnosticManager().output();

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

      #ifdef UTIL_MPI
      #ifdef MCMD_PERTURB
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
   void McSimulation::analyze(int min, int max, std::string dumpPrefix)
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

         #ifndef INTER_NOPAIR
         // Build the system CellList
         system().pairPotential().buildCellList();
         #endif

         #ifdef UTIL_DEBUG
         isValid();
         #endif

         // Initialize diagnostics (taking in molecular information).
         if (iStep_ == min) diagnosticManager().setup();

         // Sample property values
         diagnosticManager().sample(iStep_);

         // Clear out the System for the next readConfig.
         system().removeAllMolecules();

      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;

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

      iStep_ = 0;
   }

   /*
   * Get the McMove factory.
   */
   Factory<McMove>& McSimulation::mcMoveFactory()
   {  return mcMoveManagerPtr_->factory(); }

   /*
   * Check validity: return true if valid, or throw Exception.
   */
   bool McSimulation::isValid() const
   {
      Simulation::isValid();
      system().isValid();

      if (&system().simulation() != this) {
         UTIL_THROW("Invalid pointer to Simulation in System");
      }

      return true;
   }

}
#endif
