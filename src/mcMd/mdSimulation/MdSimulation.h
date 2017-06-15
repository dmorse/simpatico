#ifndef MCMD_MD_SIMULATION_H
#define MCMD_MD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdSystem.h"
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   using namespace Util;

   class MdAnalyzerManager;

   /**
   * A molecular dynamics simulation of a single MdSystem.
   *
   * \ingroup McMd_Simulation_Module
   */
   class MdSimulation : public Simulation
   {

   public:

      using ParamComposite::load;

      #ifdef UTIL_MPI
      /**
      * Constructor.
      */
      MdSimulation(MPI::Intracomm& communicator);
      #endif

      /**
      * Constructor.
      */
      MdSimulation();

      /**
      * Destructor.
      */
      virtual ~MdSimulation();

      /// \name Initialization
      //@{

      /**
      * Process command line options.
      *
      * Options:
      *
      *   -e  Enable echoing of the parameter file to a file.
      *
      *   -p  Enable use of a free energy perturbation. 
      *
      *   -r filename. Restart a simulation.
      *
      * When restarting a simulation, the required parameter "filename"
      * is the base name for the 3 required input files: filename.prm, 
      * filename.rst, and filename.cmd.
      *
      * \param argc number of arguments
      * \param argv vector of pointers to char* string arguments
      */
      void setOptions(int argc, char **argv);

      /**
      * Read parameter file.
      *
      * Returns and does nothing if in process of restarting.
      *
      * \param in parameter file stream
      */
      void readParam(std::istream &in);
   
      /**
      * Read parameters from the default parameter istream.
      *
      * Calls readParam(std::istream& ) internally, with a
      * default parameter file istream.  The default file is 
      * std::cin in serial mode (ifndef UTIL_MPI) and the 
      * file "n/param" for processor n in parallel mode 
      * (ifdef UTIL_MPI).
      */
      void readParam();

      /**
      * Read parameters from stream, without begin and end lines.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream &in);

      //@}
      /// \name Serialization and Restart Files
      //@{

      /**
      * Load parameter file section of archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Read a restart file.
      *
      * \param filename base file name for all restart files.
      */
      void load(const std::string& filename);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Write a restart file.
      *
      * \param filename base file name for all restart files.
      */
      void save(const std::string& filename);

      //@}
      /// \name Command Script
      //@{

      /**
      * Read and execute a single command from an input stream.
      *
      * Returns true if command string is recognized, false otherwise.
      *
      * \param command  command name string
      * \param in  command input stream
      */ 
      bool readCommand(std::string const & command, std::istream& in);

      /**
      * Read and execute commands from a specific input stream.
      * 
      * \param in command script input stream.
      */
      void readCommands(std::istream& in);

      /**
      * Read and execute commands from a default command file.
      *
      * This method opens a file using the commandFile name from
      * the FileMaster.
      */
      void readCommands();

      //@}
      /// \name Simulation and Analysis Operations
      //@{

      /**
      * Run an MD simulation of specified length.
      *
      * This method implements the main MD loop. The step counter iStep_ is
      * incremented until it reaches endStep. Upon exit, iStep_ = endStep.
      *
      * If isContinuation is false, the step counter iStep_ is initialized 
      * to zero, and analyzers and MdIntegrator are set to default initial 
      * states before entering the main loop. If isContinuation is true, no
      * initialization is done before entering the main loop.
      *  
      * Preconditions: Atomic positions, velocities, and forces must all 
      * be initialized before this function is entered. If isContinuation
      * is false, positions and forces must have been initialized by 
      * calling the readConfig() or load() method. The readconfig() 
      * method initializes velocities only for file formats that contain 
      * velocities. Otherwise, velocities can be initialized by the
      * system().setBoltzmannVelocities() function, which can be invoked 
      * by a THERMALIZE command in the command script.
      *
      * \param endStep         Final value of MD step counter iStep_.
      * \param isContinuation  Is this a continuation of previous run?
      */
      void simulate(int endStep, bool isContinuation = false);

      /**
      * Read and analyze a sequence of configuration files.
      *
      * This method reads and analyzes a sequence of configuration files,
      * which were normally generated by running a previous simulation using 
      * DumpConfig, and applies the sample() method of every Analyzer to
      * each such configuration. 
      *
      * The method reads files with names of the form inputPrefix() + n for 
      * integer suffixes min <= n <= max. This is consistent with the output
      * format used by DumpConfig, if the inputPrefix of the FileMaster for
      * this McAnalyzer is set to the output prefix for the dump files. 
      * 
      * \param min  integer index of first configuration file
      * \param max  integer index of last configuration file
      * \param basename  prefix string for configuration files.
      */  
      void analyzeConfigs(int min, int max, std::string basename);

      /**
      * Read and analyze a trajectory dump.
      * 
      * \param min start at this frame number
      * \param max end at this frame number
      * \param classname name of the TrajectoryReader class to use
      * \param filename  name of the trajectory file
      */
      void analyzeTrajectory(int min, int max, 
                             std::string classname, std::string filename);

      //@}
      /// \name Miscellaneous
      //@{

      /**
      * Get the MdSystem being simulated by const reference.
      */
      const MdSystem& system() const;

      /**
      * Get the MdSystem being simulated by reference.
      */
      MdSystem& system();

      /**
      * Return true if this MdSimulation valid, or throw an Exception.
      */
      virtual bool isValid() const;
 
      //@}

   private:
      
      /// System.
      MdSystem  system_;

      /// Pointer to manager for Analyzer objects.
      MdAnalyzerManager*  mdAnalyzerManagerPtr_;

      /// Restart output file name
      std::string saveFileName_;

      /// Interval for writing restart files (no output if 0)
      int saveInterval_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Is this MdSimulation in the process of restarting?
      bool isRestarting_;

   };

   // Inline method definitions

   // Get the System being simulated by const reference.
   inline const MdSystem& MdSimulation::system() const
   {  return system_; }

   // Get the System being simulated by reference.
   inline MdSystem& MdSimulation::system()
   {  return system_; }

}    
#endif
