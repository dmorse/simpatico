/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "McSimulation.h"
#include "McAnalyzerManager.h"
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/mcMoves/McMoveManager.h>
#include <mcMd/trajectory/TrajectoryReader.h>
#include <mcMd/generators/Generator.h>
#include <mcMd/generators/generatorFactory.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef UTIL_MPI
#ifdef MCMD_PERTURB
#include <mcMd/perturb/ReplicaMove.h>
#endif
#endif

#include <simp/species/Species.h>

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
   using namespace Simp;

   #ifdef UTIL_MPI
   /*
   * Constructor.
   */
   McSimulation::McSimulation(MPI::Intracomm& communicator)
    : Simulation(communicator),
      system_(),
      mcMoveManagerPtr_(0),
      mcAnalyzerManagerPtr_(0),
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
   McSimulation::McSimulation()
    : Simulation(),
      system_(),
      mcMoveManagerPtr_(0),
      mcAnalyzerManagerPtr_(0),
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

      // Create McMove and Analyzer managers
      mcMoveManagerPtr_ = new McMoveManager(*this);
      mcAnalyzerManagerPtr_ = new McAnalyzerManager(*this);

      // Pass Manager<Analyzer>* to Simulation base class.
      setAnalyzerManager(mcAnalyzerManagerPtr_);
   }

   /*
   * Destructor.
   */
   McSimulation::~McSimulation()
   {
      delete mcMoveManagerPtr_;
      delete mcAnalyzerManagerPtr_;
   }

   /*
   * Process command line options.
   */
   void McSimulation::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool rFlag = false;  // restart file
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      #ifdef MCMD_PERTURB
      bool  fflag = false;  // free energy perturbation
      #endif
      char* rarg = 0;
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f")) != -1) {
         switch (c) {
         case 'e':
            eflag = true;
            break;
         case 'r':
            rFlag = true;
            rarg  = optarg;
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            oFlag = true;
            oArg  = optarg;
            break;
         #ifdef MCMD_PERTURB
         case 'f':
           fflag = true;
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
      if (fflag) {
   
         if (rFlag) {
            std::string msg("Error: Options -r and -p are incompatible. Use -r alone. ");
            msg += "Existence of a perturbation is specified in restart file.";
            UTIL_THROW(msg.c_str());
         }
   
         // Set to expect perturbation in the param file.
         system().setExpectPerturbation();
   
         #ifdef UTIL_MPI
         Util::Log::file() << "Set to read parameters from a single file" 
                           << std::endl;
         setIoCommunicator();
         #endif
   
      }
      #endif

      // If option -p, set parameter file name
      if (pFlag) {
         if (rFlag) {
            UTIL_THROW("Cannot have both parameter and restart files");
         }
         fileMaster().setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

      // If option -r, restart
      if (rFlag) {
         //Log::file() << "Reading restart file " 
         //            << std::string(rarg) << std::endl;
         isRestarting_ = true; 
         load(std::string(rarg));
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

      // Read all species, analyzers, random number seed
      Simulation::readParameters(in);

      // Read the McSystem parameters: potential parameters, temperature etc.
      readParamComposite(in, system());

      // Read Monte Carlo Moves
      assert(mcMoveManagerPtr_);
      readParamComposite(in, *mcMoveManagerPtr_);

      // Read Analyzers (optionally)
      Analyzer::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager());

      // Parameters for writing restart files (optionally)
      saveInterval_ = 0; // default value
      readOptional<int>(in, "saveInterval", saveInterval_);
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
      loadParamComposite(ar, analyzerManager());
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
      analyzerManager().save(ar);
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
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary;
      fileMaster().openRestartIFile(filename, ar.file(), mode);
      load(ar);
      ar.file().close();

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
      Serializable::OArchive ar;
      std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary;
      fileMaster().openRestartOFile(filename, ar.file(), mode);
      save(ar);
      ar.file().close();
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
            } else {
               bool success;
               success = readCommand(command, inBuffer);
               if (!success)  {
                  Log::file() << "Error: Unknown command  " << std::endl;
                  readNext = false;
               }
            }

         }
      }
   }


   bool McSimulation::readCommand(std::string const & command, 
                                  std::istream& in)
   {
      std::string   filename;
      std::ifstream inputFile;
      std::ofstream outputFile;
      bool success = true;

      if (command == "SET_CONFIG_IO") {
         std::string classname;
         in >> classname;
         Log::file() << Str(classname, 15) << std::endl;
         system().setConfigIo(classname);
      } else
      if (command == "READ_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         fileMaster().openInputFile(filename, inputFile);
         system().readConfig(inputFile);
         inputFile.close();
      } else
      if (command == "SIMULATE") {
         int endStep;
         in >> endStep;
         Log::file() << "  " << endStep << std::endl;
         bool isContinuation = false;
         simulate(endStep, isContinuation);
      } else
      if (command == "CONTINUE") {
         if (iStep_ == 0) {
            UTIL_THROW("Attempt to continue when iStep_ == 0");
         }
         int endStep;
         in >> endStep;
         Log::file() << Int(endStep, 15) << std::endl;
         bool isContinuation = true;
         simulate(endStep, isContinuation);
      } else
      if (command == "ANALYZE_CONFIGS") {
         int min, max;
         in >> min >> max >> filename;
         Log::file() << "  " <<  min << "  " <<  max
                     << "  " <<  filename << std::endl;
         analyzeConfigs(min, max, filename);
      } else
      if (command == "ANALYZE_TRAJECTORY") {
         std::string classname;
         std::string filename;
         int min, max;
         in >> min >> max >> classname >> filename;
         Log::file() << " " << Str(classname,15) 
                     << " " << Str(filename, 15)
                     << std::endl;
         analyzeTrajectory(min, max, classname, filename);
      } else 
      if (command == "WRITE_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();
      } else
      if (command == "WRITE_PARAM") {
         in >> filename;
         Log::file() << "  " << filename << std::endl;
         fileMaster().openOutputFile(filename, outputFile);
         writeParam(outputFile);
         outputFile.close();
      } else 
      if (command == "GENERATE_MOLECULES") {
         DArray<double> diameters;
         DArray<int> capacities;
         diameters.allocate(nAtomType());
         capacities.allocate(nSpecies());

         // Parse command
         in >> system().boundary();
         Log::file() << "\n  Boundary:    " << system().boundary();
         Label capacityLabel("Capacities:");
         in >> capacityLabel;
         Log::file() << "\n  Capacities: ";
         for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
            in >> capacities[iSpecies];
            Log::file() << "  " << capacities[iSpecies];
         }
         Label diameterLabel("Diameters:");
         in >> diameterLabel;
         Log::file() << "\n  Diameters: ";
         for (int iType=0; iType < nAtomType(); iType++) {
            in >> diameters[iType];
            Log::file() << "  " << diameters[iType];
         }
         Log::file() << std::endl;

         system().generateMolecules(capacities, diameters);

      } else
      if (command == "DEFORM_CELL") {
         
         // Read in configuration from file
         in >> filename;
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
         in >> system().boundary();
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
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();

         #ifndef SIMP_NOPAIR 
         // Generate cell list
         system().pairPotential().buildCellList();
         #endif

      } else
      #ifndef UTIL_MPI
      #ifndef SIMP_NOPAIR
      if (command == "SET_PAIR") {
         std::string paramName;
         int typeId1, typeId2; 
         double value;
         in >> paramName >> typeId1 >> typeId2 >> value;
         Log::file() << "  " <<  paramName 
                     << "  " <<  typeId1 << "  " <<  typeId2
                     << "  " <<  value << std::endl;
         system().pairPotential()
                 .set(paramName, typeId1, typeId2, value);
      } else 
      #endif 
      #ifdef SIMP_BOND
      if (command == "SET_BOND") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().bondPotential().set(paramName, typeId, value);
      } else 
      #endif
      #ifdef SIMP_ANGLE
      if (command == "SET_ANGLE") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().anglePotential().set(paramName, typeId, value);
      } else 
      #endif 
      #ifdef SIMP_DIHEDRAL
      if (command == "SET_DIHEDRAL") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().dihedralPotential().set(paramName, typeId, value);
      } else 
      #endif // ifdef SIMP_DIHEDRAL
      #endif // ifndef UTIL_MPI

      // Command name not recognized
      {
         success = false;
      }
      return success;
   }

   /*
   * Read and execute commands from the default command file.
   */
   void McSimulation::readCommands()
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }


   /*
   * Run this MC simulation.
   */
   void McSimulation::simulate(int endStep, bool isContinuation)
   {
      if (!isInitialized_) {
         UTIL_THROW("McSimulation not initialized");
      }

      // Setup before main loop
      if (isContinuation) {
         Log::file() << "Restarting from iStep = " 
                     << iStep_ << std::endl;
      } else {
         iStep_ = 0;
         analyzerManager().setup();
         mcMoveManagerPtr_->setup();
      }
      int beginStep = iStep_;
      int nStep = endStep - beginStep;
      Log::file() << std::endl;
      system().positionSignal().notify();

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for ( ; iStep_ < endStep; ++iStep_) {

         // Call analyzers 
         if (Analyzer::baseInterval != 0) {
            if (iStep_ % Analyzer::baseInterval == 0) {
               if (analyzerManager().size() > 0) {
                  system().positionSignal().notify();
                  analyzerManager().sample(iStep_);
                  system().positionSignal().notify();
               }
            }
         }

         // Save restart file
         if (saveInterval_ != 0) {
            if (iStep_ % saveInterval_ == 0) {
               save(saveFileName_);
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
                  system().positionSignal().notify();
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

      // Final analyzers 
      assert(iStep_ == endStep);
      if (Analyzer::baseInterval > 0) {
         if (iStep_ % Analyzer::baseInterval == 0) {
            if (analyzerManager().size() != 0) {
               system().positionSignal().notify();
               analyzerManager().sample(iStep_);
               system().positionSignal().notify();
            }
         }
      }

      // Final save to archive
      if (saveInterval_ != 0) {
         if (iStep_ % saveInterval_ == 0) {
            save(saveFileName_);
         }
      }

      // Output results of all analyzers to output files
      if (Analyzer::baseInterval > 0) {
         analyzerManager().output();
      }

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
   void 
   McSimulation::analyzeConfigs(int min, int max, std::string basename)
   {
      // Preconditions
      UTIL_CHECK(min > 0);
      UTIL_CHECK(max > min);
      UTIL_CHECK(Analyzer::baseInterval > 0);
      UTIL_CHECK(analyzerManager().size() > 0);
      

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
         filename = basename;
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

   /*
   * Open, read and analyze a trajectory file
   */
   void McSimulation::analyzeTrajectory(int min, int max, 
                                        std::string classname, 
                                        std::string filename)
   {
      // Preconditions
      if (min < 0) UTIL_THROW("min < 0");
      if (max < 0) UTIL_THROW("max < 0");
      if (max < min) UTIL_THROW("max < min!");

      // Construct TrajectoryReader
      TrajectoryReader* trajectoryReaderPtr;
      trajectoryReaderPtr = system().trajectoryReaderFactory().factory(classname);
      if (!trajectoryReaderPtr) {
         std::string message;
         message = "Invalid TrajectoryReader class name " + classname;
         UTIL_THROW(message.c_str());
      }

      // Open trajectory file
      Log::file() << "Reading " << filename << std::endl;
      trajectoryReaderPtr->open(filename);

      // Main loop over trajectory frames
      Timer timer;
      Log::file() << "Begin main loop" << std::endl;
      bool hasFrame = true;
      timer.start();
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         hasFrame = trajectoryReaderPtr->readFrame();
         if (hasFrame) {
            #ifndef SIMP_NOPAIR
            // Build the system PairList
            system().pairPotential().buildCellList();
            #endif
            #ifdef UTIL_DEBUG
            isValid();
            #endif
            // Initialize analyzers (taking in molecular information).
            if (iStep_ == min) analyzerManager().setup();
            // Sample property values only for iStep >= min
            if (iStep_ >= min) analyzerManager().sample(iStep_);
         }
      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;
      int nFrames = iStep_ - min;

      trajectoryReaderPtr->close();
      delete trajectoryReaderPtr;

      // Output results of all analyzers to output files
      analyzerManager().output();

      // Output time 
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;
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
