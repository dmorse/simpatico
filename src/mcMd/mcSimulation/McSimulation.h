#ifndef MCMD_MC_SIMULATION_H
#define MCMD_MC_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/Simulation.h>   // base class
#include <mcMd/mcSimulation/McSystem.h>   // class member
#include <util/global.h>

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   using namespace Util;

   class McMove;
   class McMoveManager;
   class McAnalyzerManager;

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
      * \param in parameter file stream
      */
      virtual void readParam(std::istream &in);
   
      /**
      * Read parameters from the default parameter stream.
      *
      * Default parameter istream is std::cin in serial mode 
      * (ifndef UTIL_MPI) and the file "n/param" for 
      * processor n in parallel mode (ifdef UTIL_MPI).
      */
      void readParam();

      /**
      * Read parameters from a specific stream.
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

      //@}
      /// \name Simulation and Analysis Operations
      //@{

      /**
      * Run an MC simulation of specified length.
      * 
      * This method implements the main MC loop. The step counter iStep_
      * is incremented until it reaches endStep. Each step involves a
      * random selection and attempt of one Markov MC move. Upon exit,
      * iStep_ = endStep.
      *
      * If isContinuation is false, the step counter iStep_ is initialized 
      * to zero, and analyzers and mcmoves are set to default initial 
      * states before entering the main loop. If isContinuation is true, 
      * no such initialization is done for iStep_, analyzers, or the
      * MC moves.
      *  
      * \param endStep        Final value of MC step counter iStep_.
      * \param isContinuation Is this a continuation of a previous run?
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
      * format used by the DumpConfig class.
      *
      * In serial mode, the inputPrefix should be given as a path relative
      * to the directory in which the program is executed. The inputPrefix 
      * of the simulation FileMaster is not prepended to the dump Prefix.  
      *
      * In parallel mode, for processor with MPI rank m, the path "m/" is 
      * prepended to the basename, so that all files associated with this
      * processor are in this directory, but no inputPrefix is added after
      * the string "m/".
      *
      * \param min  integer suffix of first configuration file name
      * \param max  integer suffix of last configuration file name
      * \param basename  root name for dump files (without integer suffix)
      */  
      void analyzeConfigs(int min, int max, std::string basename);

      /**
      * Read and analyze a trajectory file.
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
      McSystem       system_;
   
      /// Pointer to Manager for Monte Carlo moves.
      McMoveManager* mcMoveManagerPtr_;

      /// Pointer to Manager for analyzers.
      McAnalyzerManager* mcAnalyzerManagerPtr_;

      /// Pointer to parameter file passed to readParameters(istream&)
      std::istream*   paramFilePtr_;

      /// Restart output file name
      std::string saveFileName_;

      /// Interval for writing restart files (no output if 0)
      int saveInterval_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Is this McSimulation in the process of restarting?
      bool isRestarting_;

   }; 

   // Inline Methods

   /* 
   * Get the McSystem.
   */
   inline McSystem& McSimulation::system()
   { return system_; }

   /* 
   * Get a const ref to the McSystem.
   */
   inline const McSystem& McSimulation::system() const
   { return system_; }

   /* 
   * Get the McMoveManager (protected).
   */
   inline McMoveManager& McSimulation::mcMoveManager()
   {
      assert(mcMoveManagerPtr_);  
      return *mcMoveManagerPtr_; 
   }

}    
#endif
