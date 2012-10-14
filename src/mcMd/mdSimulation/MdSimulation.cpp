#ifndef MCMD_MD_SIMULATION_CPP
#define MCMD_MD_SIMULATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdDiagnosticManager.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/trajectoryIos/TrajectoryIo.h>
#include <mcMd/species/Species.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/misc/Log.h>
#include <util/random/serialize.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <util/misc/Timer.h>

#include <sstream>
#include <iostream>
#include <unistd.h>

namespace McMd
{

   using namespace Util;

   #ifdef UTIL_MPI
   /* 
   * Constructor.
   */
   MdSimulation::MdSimulation(MPI::Intracomm& communicator)
    : Simulation(communicator),
      system_(),
      mdDiagnosticManagerPtr_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("MdSimulation"); 
      system_.setId(0);
      system_.setSimulation(*this);
      system_.setFileMaster(fileMaster());
      mdDiagnosticManagerPtr_ = new MdDiagnosticManager(*this);
      assert(mdDiagnosticManagerPtr_);
      setDiagnosticManager(mdDiagnosticManagerPtr_);
   }
   #endif

   /* 
   * Constructor.
   */
   MdSimulation::MdSimulation()
    : Simulation(),
      system_(),
      mdDiagnosticManagerPtr_(0),
      isInitialized_(false),
      isRestarting_(false)
   {
      setClassName("MdSimulation"); 
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
   * Process command line options.
   */
   void MdSimulation::setOptions(int argc, char **argv)
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
      if (pflag) {
   
         // Set to expect perturbation in the param file.
         system().setExpectPerturbation();
   
         #ifdef UTIL_MPI
         Util::Log::file() << "Set to read parameters from a single file" 
                           << std::endl;
         setParamCommunicator();
         #endif
   
      }
      #endif
      if (rflag) {
         std::cout << "Reading restart" << std::endl;
         std::cout << "Base file name " << std::string(rarg) << std::endl;
         isRestarting_ = true; 
         readRestart(std::string(rarg));
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
      readParamComposite(in, diagnosticManager());

      isValid();
      isInitialized_ = true;
   }
 
   /*
   * Read default parameter file.
   */
   void MdSimulation::readParam()
   {  readParam(fileMaster().paramFile()); }

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

            if (command == "FINISH") {
               Log::file() << std::endl;
               readNext = false;
            } else 
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
               // Generate pair list
               system().pairPotential().buildPairList();
               #endif

            } else
            #ifndef UTIL_MPI
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
            #endif
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
            #endif
            #endif
            {
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
         diagnosticManager().setup();
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

         #ifdef INTER_NOPAIR
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

      #ifndef INTER_NOPAIR
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

         #ifndef INTER_NOPAIR
         // Build the system CellList
         system().pairPotential().buildPairList();
         #endif

         #ifndef NDEBUG
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

         #ifndef INTER_NOPAIR
         // Build the system CellList
         system().pairPotential().buildPairList();
         #endif

         #ifndef NDEBUG
         isValid();
         #endif

         // Initialize diagnostics (taking in molecular information).
         if (iStep_ == min) diagnosticManager().setup();

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
      // readRestart
      if (isInitialized_) {
         UTIL_THROW("Error: Called readRestart when already initialized");
      }
      if (!isRestarting_) {
         UTIL_THROW("Error: Called readRestart without restart option");
      }

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

      // Set command (*.cmd) file
      std::string commandFileName = filename + ".cmd";
      fileMaster().setCommandFileName(commandFileName);

      isInitialized_ = true;
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
