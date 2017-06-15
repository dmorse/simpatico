/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdAnalyzerManager.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/generators/Generator.h>
#include <mcMd/generators/generatorFactory.h>
#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/trajectory/TrajectoryReader.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#include <simp/species/Species.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/misc/Log.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <util/misc/Timer.h>

#include <sstream>
#include <iostream>
#include <unistd.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   #ifdef UTIL_MPI
   /* 
   * Constructor.
   */
   MdSimulation::MdSimulation(MPI::Intracomm& communicator)
    : Simulation(communicator),
      system_(),
      mdAnalyzerManagerPtr_(0),
      saveFileName_(),
      saveInterval_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("MdSimulation"); 
      system_.setId(0);
      system_.setSimulation(*this);
      system_.setFileMaster(fileMaster());
      mdAnalyzerManagerPtr_ = new MdAnalyzerManager(*this);
      assert(mdAnalyzerManagerPtr_);
      setAnalyzerManager(mdAnalyzerManagerPtr_);
   }
   #endif

   /* 
   * Constructor.
   */
   MdSimulation::MdSimulation()
    : Simulation(),
      system_(),
      mdAnalyzerManagerPtr_(0),
      saveFileName_(),
      saveInterval_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("MdSimulation"); 
      system_.setId(0);
      system_.setSimulation(*this);
      system_.setFileMaster(fileMaster());
      mdAnalyzerManagerPtr_ = new MdAnalyzerManager(*this);
      assert(mdAnalyzerManagerPtr_);
      setAnalyzerManager(mdAnalyzerManagerPtr_);
   }

   /* 
   * Destructor.
   */
   MdSimulation::~MdSimulation()
   {
      assert(mdAnalyzerManagerPtr_);
      if (mdAnalyzerManagerPtr_) {
         delete mdAnalyzerManagerPtr_;
         mdAnalyzerManagerPtr_ = 0;
      }
   }

   /*
   * Process command line options.
   */
   void MdSimulation::setOptions(int argc, char **argv)
   {
      bool  eflag = false;  // echo
      bool  rFlag = false;  // restart file
      bool  pFlag = false;  // param file 
      bool  cFlag = false;  // command file 
      bool  iFlag = false;  // input prefix
      bool  oFlag = false;  // output prefix
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
           std::cout << "Unknown option -" << optopt << std::endl;
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
            std::string msg("Error: Options -r and -f are incompatible.");
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
         // Log::file() << "Reading restart file "
         //             << std::string(rarg) << std::endl;
         isRestarting_ = true; 
         load(std::string(rarg));
      }

   }

   /* 
   * Read parameters from file.
   */
   void MdSimulation::readParameters(std::istream &in)
   { 
      if (isInitialized_) {
         UTIL_THROW("Error: Called readParam when already initialized");
      }

      Simulation::readParameters(in); 

      readParamComposite(in, system_);
      Analyzer::baseInterval = 0; 
      readParamCompositeOptional(in, analyzerManager());

      // Parameters for writing restart files
      saveInterval_ = 0;
      readOptional<int>(in, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         read<std::string>(in, "saveFileName", saveFileName_);
      }

      isValid();
      isInitialized_ = true;
   }
 
   /*
   * Read default parameter file.
   */
   void MdSimulation::readParam()
   {
      readParam(fileMaster().paramFile()); 
   }

   /*
   * Read parameter block, including begin and end.
   */
   void MdSimulation::readParam(std::istream &in)
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
   * Load parameters from an archive.
   */
   void MdSimulation::loadParameters(Serializable::IArchive& ar)
   { 
      if (isInitialized_) {
         UTIL_THROW("Error: Called loadParameters when already initialized");
      }

      Simulation::loadParameters(ar); 
      loadParamComposite(ar, system_); 
      loadParamComposite(ar, analyzerManager());
      loadParameter<int>(ar, "saveInterval", saveInterval_);
      if (saveInterval_ > 0) {
         loadParameter<std::string>(ar, "saveFileName", saveFileName_);
      }

      system_.loadConfig(ar);
      ar & iStep_;

      isValid();
      isInitialized_ = true;
   }
 
   /* 
   * Save parameters to an archive.
   */
   void MdSimulation::save(Serializable::OArchive& ar)
   { 
      Simulation::save(ar); 
      system_.saveParameters(ar);
      analyzerManager().save(ar);
      ar << saveInterval_;
      if (saveInterval_ > 0) {
         ar << saveFileName_;
      }
      system_.saveConfig(ar);
      ar & iStep_;
   }
 
   /*
   * Read and implement commands in an input script.
   */
   void MdSimulation::readCommands(std::istream &in)
   {
      std::string   command;
      #ifndef UTIL_MPI
      std::istream&     inBuffer = in;
      #else
      std::stringstream inBuffer;
      std::string       line;
      #endif

      bool readNext = true;
      while (readNext) {

         #ifdef UTIL_MPI
         if (!hasIoCommunicator() || isIoProcessor()) {
            getNextLine(in, line);
         } 
         if (hasIoCommunicator()) {
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

   /*
   * Read and implement commands from the default command file.
   */
   void MdSimulation::readCommands()
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }

   bool MdSimulation::readCommand(std::string const & command, 
                                  std::istream& in)
   {
      std::string   filename;
      std::ifstream inputFile;
      std::ofstream outputFile;
      bool success = true;

      if (command == "READ_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         fileMaster().openInputFile(filename, inputFile);
         std::cout << "Opened config file" << std::endl;
         system().readConfig(inputFile);
         std::cout << "Finished reading config file" << std::endl;
         inputFile.close();
      } else
      if (command == "THERMALIZE") {
         double temperature;
         in >> temperature;
         Log::file() << Dbl(temperature, 15, 6) << std::endl;
         system().setBoltzmannVelocities(temperature);
         system().removeDriftVelocity();
      } else
      if (command == "SIMULATE") {
         int endStep;
         in >> endStep;
         Log::file() << Int(endStep, 15) << std::endl;
         bool isContinuation = false;
         simulate(endStep, isContinuation);
      } else
      if (command == "CONTINUE") {
         if (iStep_ == 0) {
            UTIL_THROW("Attempt to continue simulation when iStep_ == 0");
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
         Log::file() <<  Int(min, 15) << Int(max, 15)
                     <<  Str(filename, 20) << std::endl;
         analyzeConfigs(min, max, filename);
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
         Log::file() << Str(filename, 15) << std::endl;
         fileMaster().openOutputFile(filename, outputFile);
         writeParam(outputFile);
         outputFile.close();
      } else
      if (command == "SET_CONFIG_IO") {
         std::string classname;
         in >> classname;
         Log::file() << Str(classname, 15) << std::endl;
         system().setConfigIo(classname);
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
      if (command == "GENERATE_MOLECULES") {
         DArray<double> diameters;
         DArray<int> capacities;
         diameters.allocate(nAtomType());
         capacities.allocate(nSpecies());

         // Parse command
         in >> system().boundary();
         Log::file() << "  " << system().boundary();
         Label capacityLabel("Capacities:");
         in >> capacityLabel;
         for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
            in >> capacities[iSpecies];
            Log::file() << "  " << capacities[iSpecies];
         }
         Label diameterLabel("Diameters:");
         in >> diameterLabel;
         for (int iType=0; iType < nAtomType(); iType++) {
            in >> diameters[iType];
            Log::file() << "  " << diameters[iType];
         }
         Log::file() << std::endl;

         system().generateMolecules(capacities, diameters);
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
      #endif 
      #endif
      {
         success = false;
      }
      return success;
   }

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
         system().shiftAtoms();
         #ifdef UTIL_COULOMB
         if (system().hasCoulombPotential()) {
            system().coulombPotential().makeWaves();
         }
         #endif
         system().calculateForces();
         analyzerManager().setup();
         system_.mdIntegrator().setup();
      }
      int beginStep = iStep_;
      int nStep = endStep - beginStep;

      #ifdef SIMP_NOPAIR
      // When the pair potential is disabled, require that
      // Analyzer::baseInterval > 0 to guarantee periodic 
      // shifting of atomic positions into primary cell.
      UTIL_CHECK(Analyzer::baseInterval > 0);
      #endif

      // Main loop 
      Log::file() << std::endl;
      timer.start();
      for ( ; iStep_ < endStep; ++iStep_) {

         // Sample scheduled analyzers
         if (Analyzer::baseInterval > 0) {
            if (iStep_ % Analyzer::baseInterval == 0) {
               system().shiftAtoms();
               analyzerManager().sample(iStep_);
            }
         }

         if (saveInterval_ > 0) {
            if (iStep_ % saveInterval_ == 0) {
               save(saveFileName_);
            }
         }

         // Take one MD step with the MdIntegrator
         system_.mdIntegrator().step();
      }
      timer.stop();
      double time  = timer.time();
      double rstep = double(nStep);

      // Shift final atomic positions 
      system().shiftAtoms();

      // Final analyzers and restart
      if (Analyzer::baseInterval > 0) {
         if (iStep_ % Analyzer::baseInterval == 0) {
            analyzerManager().sample(iStep_);
         }
      }
      if (saveInterval_ > 0) {
         if (iStep_ % saveInterval_ == 0) {
            save(saveFileName_);
         }
      }

      // Final analyzer output
      analyzerManager().output();

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

      #ifndef SIMP_NOPAIR
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
   void 
   MdSimulation::analyzeConfigs(int min, int max, std::string basename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max > min);
      UTIL_CHECK(Analyzer::baseInterval);
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
         system().pairPotential().buildPairList();
         #endif

         #ifndef NDEBUG
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
      Log::file()<< "end main loop" << std::endl;

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

   }

   /*
   * Open, read and analyze a trajectory file
   */
   void MdSimulation::analyzeTrajectory(int min, int max, 
                                        std::string classname, 
                                        std::string filename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max > min);
      UTIL_CHECK(Analyzer::baseInterval);
      UTIL_CHECK(analyzerManager().size() > 0);

      // Construct TrajectoryReader
      TrajectoryReader* trajectoryReaderPtr;
      trajectoryReaderPtr 
             = system().trajectoryReaderFactory().factory(classname);
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
            system().pairPotential().buildPairList();
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

   /**
   * Open, load, and close binary restart file. 
   */
   void MdSimulation::load(const std::string& filename)
   {
      // Precondition
      if (isInitialized_) {
         UTIL_THROW("Error: Called load when already initialized");
      }

      // Load state from archive
      Serializable::IArchive ar;
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary;
      fileMaster().openRestartIFile(filename, ar.file(), mode);
      load(ar);
      ar.file().close();

      isInitialized_ = true;
      isRestarting_  = true;
   }

   /*
   * Write a restart file with specified filename, if at interval.
   */
   void MdSimulation::save(const std::string& filename)
   {
      Serializable::OArchive ar;
      std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary;
      fileMaster().openRestartOFile(filename, ar.file(), mode);
      save(ar);
      ar.file().close();
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
