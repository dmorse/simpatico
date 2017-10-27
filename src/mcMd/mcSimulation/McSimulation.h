#ifndef MCMD_MC_SIMULATION_H
#define MCMD_MC_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/Simulation.h>           // base class
#include <mcMd/mcSimulation/McSystem.h>           // member
#include <mcMd/mcSimulation/McAnalyzerManager.h>  // member
#include <mcMd/mcSimulation/McCommandManager.h>   // member
#include <mcMd/mcMoves/McMoveManager.h>           // member
#include <util/global.h>

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   using namespace Util;

   class McMove;

   /**
   * A Monte-Carlo simulation of one McSystem.
   *
   * \ingroup McMd_Simulation_Module
   */
   class McSimulation : public Simulation
   {

   public:

      using ParamComposite::load;

      #ifdef UTIL_MPI
      /**
      * Constructor.
      */
      McSimulation(MPI::Intracomm& communicator);
      #endif

      /**
      * Constructor.
      */
      McSimulation();

      /**
      * Destructor.
      */
      virtual ~McSimulation();

      /// \name Initialization
      //@{

      /**
      * Process command line options.
      *
      * Main options:
      *  
      *   -q           Query: Print enabled/disabled features
      *
      *   -e           Echo: Enable echoing of parameter file
      *
      *   -r filename  Restart: restart from specified file
      *
      *   -p filename  Parameter: Specify a parameter file
      *
      *   -c filename  Command: Specify a command file
      * 
      *   -i path      Input: Specify path prefix for input files
      *
      *   -o path      Input: Specify path prefix for output files
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
      * file named n/<filename>, where <filename> is either the
      * argument of the -p command line option, if invoked with 
      * that option, or the default string <filename> = "param".
      *
      * \pre: Call after setOptions().
      */
      void readParam();

      /**
      * Read specified parameter file.
      *
      * Returns and does nothing if in process of restarting
      * (i.e., if the main program was invoked with -r option).
      * This calls readParameters(std::istream& ) internally.
      *
      * \pre: Call after setOptions().
      *
      * \param in parameter file stream
      */
      virtual void readParam(std::istream &in);
   
      /**
      * Read body of parameter block from a specific file.
      *
      * \param in parameter file input stream.
      */
      virtual void readParameters(std::istream &in);

      //@}
      /// \name Serialization and Restarting
      //@{

      /**
      * Load parameters from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Read restart file to continue simulation.
      *
      * \param filename name of input restart file
      */
      void load(const std::string& filename);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Write internal state to a restart file.
      *
      * \param filename name of output restart file
      */
      void save(const std::string& filename);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, unsigned int version);

      //@}
      /// \name Command Script Interface
      //@{

      /**
      * Read and execute commands from a command file input stream.
      * 
      * \param in  command file input stream
      */
      void readCommands(std::istream& in);

      /**
      * Read and execute commands from the default parameter file.
      * 
      * This method opens the file with the name commandFile read
      * by the FileMaster.
      */
      void readCommands();

      /**
      * Read and execute a single command from an input stream.
      *
      * Usage: The capitalized command name string must have been 
      * read from istream "in" and passed as the "command" argument. 
      * If the command name is recognized, any additional arguments 
      * are read from stream "in", the command is executed, and a
      * value of true is returned. A value of false is returned
      * iff the command name string is not recognized.
      * 
      * Implementation: Calls CommandManager::readCommand().
      *
      * \param command  command name string
      * \param in  command input stream
      */ 
      bool readCommand(std::string command, std::istream& in);

      //@}
      /// \name Simulation and Analysis Operations
      //@{

      /**
      * Run an MC simulation of specified length.
      * 
      * This method implements the main MC loop. The step counter iStep
      * is incremented until it reaches endStep. Each step involves a
      * random selection and attempt of one Markov MC move. Upon exit,
      * iStep_ = endStep.
      *
      * If isContinuation is false, the step counter iStep is initialized 
      * to zero, and analyzers and mcmoves are set to default initial 
      * states before entering the main loop. If isContinuation is true, 
      * no such initialization is done for iStep_, analyzers, or the
      * MC moves.
      *  
      * \param endStep  Final value of MC step counter iStep_.
      * \param isContinuation  Is this a continuation of a previous run?
      */
      void simulate(int endStep, bool isContinuation = false);

      /**
      * Read and analyze a sequence of configuration files.
      *
      * This method reads and analyzes a sequence of configuration files,
      * which were normally generated by running a previous simulation 
      * using ConfigWriter, and applies the sample() method of every 
      * Analyzer to each such configuration. 
      *
      * The method reads a sequence of configuration files with names of 
      * the form inputPrefix + basename + n for integer suffixes in the
      * range min <= n <= max. This is consistent with the output format
      * format used by the WriteConfig class. The inputPrefix used in 
      * an analysis simulation is often a directory name, with a trailing
      * directory separator "/", that is equal to the outputPrefix used
      * in the earlier simulation run.
      *
      * In parallel mode, for processor with MPI rank m, the path "m/" 
      * is prepended to the fileMaster input prefix, so that paths to 
      * all files associated with processor m begin with the string
      * "m/inputPrefix" + basename.
      *
      * \param min  integer suffix of first configuration file name
      * \param max  integer suffix of last configuration file name
      * \param basename  root name for dump files (without integer suffix)
      */  
      void analyzeConfigs(int min, int max, std::string basename);

      /**
      * Read and analyze a trajectory file.
      * 
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory 
      * file with a path of the form inputPrefix + filename. 
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader class to use
      * \param filename  name of the trajectory file
      */
      void analyzeTrajectory(int min, int max, 
                             std::string classname, std::string filename);

      //@}
      /// \name Miscellaneous
      //@{

      /**
      * Get the McSystem by reference.
      */
      McSystem& system();

      /**
      * Get the McSystem by const refererence.
      */
      const McSystem& system() const;

      /**
      * Get the McMove factory by reference.
      */
      Factory<McMove>& mcMoveFactory();

      /**
      * Return true if valid, or throw an Exception. 
      */
      virtual bool isValid() const;

      //@}

   protected:

      /**
      * Get the McMoveManager by reference.
      */
      McMoveManager& mcMoveManager();

   private:
   
      /// System.
      McSystem system_;
   
      /// Manager for Monte Carlo moves.
      McMoveManager mcMoveManager_;

      /// Manager for Analyzer objects.
      McAnalyzerManager mcAnalyzerManager_;

      /// Manager for Command objects.
      McCommandManager mcCommandManager_;

      /// Pointer to parameter file passed to readParameters(istream&)
      std::istream* paramFilePtr_;

      /// Restart output file name
      std::string saveFileName_;

      /// Interval for writing restart files (no output if 0)
      int saveInterval_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Is this McSimulation in the process of restarting?
      bool isRestarting_;

   }; 

   // Inline member functions

   /* 
   * Get the McSystem by reference.
   */
   inline McSystem& McSimulation::system()
   {  return system_; }

   /* 
   * Get the McSystem by const reference.
   */
   inline const McSystem& McSimulation::system() const
   {  return system_; }

   /* 
   * Get the McMoveManager (protected).
   */
   inline 
   McMoveManager& McSimulation::mcMoveManager()
   {  return mcMoveManager_; }

}    
#endif
