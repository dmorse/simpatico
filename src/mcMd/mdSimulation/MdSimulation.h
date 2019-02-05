#ifndef MCMD_MD_SIMULATION_H
#define MCMD_MD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/Simulation.h>     // base class
#include "MdSystem.h"                       // member
#include "MdAnalyzerManager.h"              // member
#include "MdCommandManager.h"               // member

namespace McMd
{

   using namespace Util;

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
      MdSimulation(MPI_Comm communicator);
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
      * Main options:
      *  
      *   -q   Query: Print list of enabled/disabled features
      *
      *   -e   Echo: Enable echoing of parameter file on read
      *
      *   -r filename.  Restart: restart from specified file
      *
      *   -p filename.  Parameter: Specify a parameter file
      *
      *   -c filename.  Command: Specify a command file
      * 
      *   -i filename.  Input: Specify path prefix for input files
      *
      *   -o filename.  Input: Specify path prefix for output files
      *
      * The -p and -r options are mutually exclusive: When a
      * simulation is restarted, all information required from
      * a parameter file is in the restart file. 
      *
      * \param argc number of arguments
      * \param argv vector of pointers to char* string arguments
      */
      void setOptions(int argc, char **argv);

      /**
      * Read parameters from the default parameter istream.
      *
      * Calls readParam(std::istream&) internally, with a
      * default parameter file istream given by the return
      * value of FileMaster::paramFile(). 
      *
      * Single parameter file: If compiled as a serial program
      * (ifndef UTIL_MPI) or as a parallel program in mode that
      * uses a single parameter file (i.e., with option -f), 
      * the parameter file name is the argument passed to the 
      * -p command line option, if the main program is invoked 
      * with the -p option, or the parameter file is read from
      * standard input, std::cin, if not invoked with a -p
      * option.
      * 
      * Multiple parameter files: If compiled as a parallel
      * program (ifdef UITL_MPI) and used in a mode with 
      * separate parameter files for independent simulations
      * (i.e., without the -f option), then the parameter file
      * for the simulation performed by processor number n is
      * file named n/filename, where "filename" is either the
      * argument of the -p command line option, if invoked with 
      * that option, or the default string filename = "param".
      *
      * \pre: Call after setOptions().
      */
      void readParam();

      /**
      * Read parameter file.
      *
      * Returns and does nothing if in process of restarting
      * (i.e., if the main program was invoked with -r option).
      *
      * \pre: Call after setOptions().
      *
      * \param in parameter file stream
      */
      void readParam(std::istream &in);
   
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
      * Read a restart file.
      *
      * \param filename base file name for all restart files.
      */
      void load(const std::string& filename);

      /**
      * Load parameter file section of archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Write a restart file.
      *
      * Calls save(Serializable::OArchive) internally. 
      *
      * \param filename restart file name
      */
      void save(const std::string& filename);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      //@}
      /// \name Command Script
      //@{

      /**
      * Read and execute commands from a default command file.
      *
      * This function calls readCommands(fileMaster().commandFile()).
      */
      void readCommands();

      /**
      * Read and execute commands from a specific input stream.
      * 
      * \param in command script input stream.
      */
      void readCommands(std::istream& in);

      /**
      * Read and execute a single command from an input stream.
      *
      * Usage: The capitalized command name must have been read
      * from istream "in" and passed as the "command" argument. 
      * If the command name is recognized, any additional arguments 
      * are read from stream "in", the command is executed, and
      * a value of true is returned. A value of false is returned
      * iff the command name string is not recognized.
      * 
      * Calls commandManager().readCommand(command, in) internally.
      *
      * \param command  command name string
      * \param in  command input stream
      */ 
      bool readCommand(std::string command, std::istream& in);

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
      * This function reads and analyzes a sequence of configuration files
      * that were normally generated by running a previous simulation using 
      * the ConfigWriter analyzer to periodically dump configurations. This
      * function applies the sample() method of every Analyzer to each such
      * stored configuration. 
      * 
      * The function reads a sequence of configuration files with names of 
      * the form inputPrefix + basename + n , for integer suffixes in the
      * range min <= n <= max. This is consistent with the output format
      * format used by the WriteConfig class. The inputPrefix used in 
      * an analysis simulation is often a directory name, with a trailing
      * directory separator "/", that is the same as the outputPrefix used
      * in the earlier simulation run.
      *
      * In parallel mode, for processor with MPI rank m, the path "m/" 
      * is prepended to the fileMaster input prefix, so that paths to 
      * all files associated with processor m begin with the string
      * "m/inputPrefix" + basename.
      * 
      * \param min  integer index of first configuration file
      * \param max  integer index of last configuration file
      * \param basename  prefix string for configuration files.
      */  
      void analyzeConfigs(int min, int max, std::string basename);

      /**
      * Read and analyze a trajectory file.
      * 
      * This function reads and analyzes a trajectory file, which is a
      * single file containing a sequence of configuration snapshots. 
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
      MdSystem system_;

      /// Manager for Analyzer objects.
      MdAnalyzerManager mdAnalyzerManager_;

      /// Manager for Command objects.
      MdCommandManager mdCommandManager_;

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
