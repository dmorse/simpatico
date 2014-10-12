#ifndef UTIL_FILE_MASTER_H
#define UTIL_FILE_MASTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class

#include <ios>
#include <fstream>
#include <string>

namespace Util
{

   /**
   * A FileMaster manages input and output files for a simulation.
   *
   * A FileMaster manages a set of input and output files for a  
   * molecular simulation. It provides methods to open four different
   * types of file, which are placed in different locations within a 
   * standard directory structure. These file types are:
   *
   *   - Control files: parameter and command input files
   *
   *   - Restart files, which can be opened for input or output
   *
   *   - Input data files (e.g., input configuration files)
   *
   *   - Output data files (e.g., trajectories or analysis outputs)
   *
   * Somewhat different conventions are used for file paths in
   * simulations of a single system and for parallel simulations 
   * of multiple systems, as discussed below. 
   * 
   * Simulations of a single system
   * ------------------------------
   *
   * In simulations of a single system, paths to the different types
   * of file are constructed by prepending some combination of the
   * following elements before a base file name:
   *
   *  - Root directory prefix: This is the path from the current
   *    working directory to the root level directory for all files
   *    associated with a simulation run. It defaults to an empty
   *    string. It may be set by calling the setRootPrefix() function.
   *
   *  - Input prefix: This string is prepended to the base name of
   *    all input data files. 
   * 
   *  - Output prefix: This string is prepended to the base name of
   *    all output data files. 
   * 
   * The root directory prefix must be a path to a directory and
   * thus must either be an empty string or must end with the 
   * directory separator "/". The input and output prefix are 
   * often also chosen to be directory paths, in order to place 
   * input and output files in different subdirectories (e.g., 
   * "in/" and "out/"), but can also be strings that are 
   * prepended to the base name for all input and output files.
   *
   * In simulations of a single systems, paths to both control
   * and restart files are constructed by the openControFile(),
   * openRestartIFile() and openRestartOFile() by concatenating
   * the root prefix (if any) to a base file name that is passed 
   * to the appropriate function. Paths to input and output files
   * are constructed by the openInputFile() and openOutputFile()
   * functions by concatenating the root prefix, if any, to the 
   * input or output prefix, as appropriate, and a base file 
   * name that is passed as as argument. 
   *
   * Simulations of a multiple systems
   * ---------------------------------
   * A parallel simulation of multiple systems use a directory
   * structure in which the root directory (specified by the root
   * directory prefix) contains a set of numbered subdirectories
   * named "0/", "1/", "2/", etc.. Each such numbered directory
   * contains input, output and restart files that are specific 
   * to a particular system, and will be referred to in what 
   * follows as the system directory for that system. The integer
   * index for the system associated with a specific MPI processor
   * must be set by the setDirectoryId(int id) method.
   *
   * In all simulations of multiple systems, paths to input and
   * input and output data files are constructed by the functions
   * openInputFile() and openOutputFile() by concatenating the
   * system directory path, an input or output prefix, and a base
   * file name. Paths to restart files for each system are 
   * constructed by concatentaing the system directory path and
   * a base file name, thus placing these files in the system
   * directory. 
   *
   * Two different modes of operation are possible for simulations
   * of multiple systems, which differ in the treatment of control
   * files. In "independent" mode, simulations of multiple systems 
   * are assumed to be completely independent and to require separate
   * parameter and command control files for each system. In this
   * case, the openControlFile() system opens files in the numbered
   * system directory. In "replicated" model, all simulations are
   * controlled by a single parameter file and a single control
   * file, which are opened in the root directory. "Independent"
   * mode is enabled by default, and "replicated" mode may be chosen
   * by calling the function "setCommonControl()" before opening any
   * any control files. 
   *   
   * Functions for Parameter and Command Files
   * -----------------------------------------
   * 
   * Base file names to the parameter and command files can be set by 
   * calling the setParamFileName() or setCommandFileName() functions,
   * respectively. If the setCommandFileName() function is not called
   * before readParameters(), the readParameters() function will
   * attempt to read a "commandFileName" parameter from the FileMaster
   * parameter file block. After parameter and command file names are
   * set, references to the parameter and command files may be 
   * retrieved by calling the paramFile() and commandFile() functions. 
   * Each of these functions calls openControlFile() function internally 
   * the first time it is invoked.
   *
   * \ingroup Misc_Module
   */
   class FileMaster : public ParamComposite
   {
   
   public:

      /**
      * Constructor.
      */
      FileMaster();
 
      /**
      * Copy constructor.
      *
      * \param copy FileMaster object to be copied
      */
      FileMaster(const FileMaster& copy);
 
      /**
      * Destructor.
      */
      virtual ~FileMaster();

      /// \name Initialization
      //@{

      /**
      * Set the path from current directory to root directory.
      *
      * \param rootPrefix root prefix for all file names.
      */
      void setRootPrefix(const std::string& rootPrefix);

      /**
      * Set an integer directory identifier for this processor.
      *
      * This method should be called only for simulations of multiple
      * systems, to set the integer index of the physical system
      * associated with this processor.  After calling 
      * setDirectoryId(n) with an integer n, a prefix "n/" will be 
      * prepended to the paths of input and output files associated
      * with that system. 
      *
      * \param directoryId integer subdirectory name
      */
      void setDirectoryId(int directoryId);

      /**
      * Set to use single param and command file for control.
      */ 
      void setCommonControl();

      /**
      * Set the parameter file name.
      *
      * \param paramFileName name of parameter file
      */
      void setParamFileName(const std::string& paramFileName);

      /**
      * Set the command file name.
      *
      * \param commandFileName name of command file
      */
      void setCommandFileName(const std::string& commandFileName);

      /**
      * Set the input file prefix string.
      *
      * \param inputPrefix input file prefix string
      */
      void setInputPrefix(const std::string& inputPrefix);
 
      /**
      * Set the output file prefix string.
      *
      * \param outputPrefix output file prefix string
      */
      void setOutputPrefix(const std::string& outputPrefix);
 
      /**
      * Read parameter file.
      *
      * Reads the inputPrefix and outputPrefix string variables.
      *
      * \param in pararameter file input stream
      */
      virtual void readParameters(std::istream& in);

      //@}
      /// \name Serialization
      //@{

      /**
      * Load internal state from file.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to file. 
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
      //@}
      /// \name File Opening
      //@{
      
      /**
      * Open an input file with a known path and open mode.
      *
      * Adds error checking to C++ ifstream::open function.
      *
      * \param  name  complete file path
      * \param  in    ifstream object to associated with a file
      * \param  mode  read mode 
      */
      void open(const std::string& name, std::ifstream& in, 
                std::ios_base::openmode mode = std::ios_base::in) const;

      /**
      * Open an output file with a known path and open mode.
      *
      * Add error checking to C++ ofstream::open function.
      *
      * \param  name  complete file path
      * \param  out   ofstream object to associated with a file
      * \param  mode  write mode
      */
      void open(const std::string& name, std::ofstream& in, 
                std::ios_base::openmode mode = std::ios_base::out) const;

      /**
      * Open an input parameter or command file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id 
      * has not been set and the setCommonControl() function has 
      * not been called.
      *
      * \param  name  base file name, without any prefix
      * \param  in    ifstream object to open
      */
      void openControlFile(const std::string& name, 
                           std::ifstream& in) const;

      /**
      * Open an input restart dump file for reading.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * \param name  base file name, without any prefix or extension
      * \param in  ifstream object to open
      * \param mode open mode
      */
      void 
      openRestartIFile(const std::string& name, std::ifstream& in,
                       std::ios_base::openmode mode = std::ios_base::in) 
      const;

      /**
      * Open an output restart file for writing.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name 
      *
      * \param name base file name
      * \param out ofstream object to open
      * \param mode open mode
      */
      void 
      openRestartOFile(const std::string& name, std::ofstream& out,
                       std::ios_base::openmode mode = std::ios_base::out) 
      const;

      /**
      * Open an input file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + inputPrefix + filename.
      *
      * \param  filename  file name, without any prefix
      * \param  in  ifstream object to associated with a file
      * \param  mode  bit mask that specifies opening mode
      */
      void 
      openInputFile(const std::string& filename, std::ifstream& in, 
                    std::ios_base::openmode mode = std::ios_base::in) 
      const;

      /**
      * Open an output file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + outputPrefix + filename.
      *
      * \param  filename  file name, without any prefix
      * \param  out  ofstream object to associated with a file
      * \param  mode  bit mask that specifies opening mode
      */
      void 
      openOutputFile(const std::string& filename, std::ofstream& out, 
                     std::ios_base::openmode mode = std::ios_base::out) 
      const;

      //@}
      /// \name Misc Accessors
      //@{
      
      /**
      * Is set for common param and command files?
      */
      bool isCommonControl() const;

      /**
      * Return the param file name, if any.
      */
      std::string paramFileName() const;

      /**
      * Return the command file name.
      *
      * The base name of the command file is read from the parameter file
      * by the readParameters() method. 
      */
      std::string commandFileName() const;

      /**
      * Get a default parameter stream by reference.
      *
      * If setDirectoryId() has not been called, of if setCommonControl()
      * has been called, this method returns std::cin.
      * 
      * If setDirectoryId() has been called and setCommonControld() has
      * not, this method returns a reference to a file "n/param". This 
      * file is opened for reading the first time it is returned by this 
      * function.
      */
      std::istream& paramFile();

      /**
      * Get the command input stream by reference.
      *
      * If the commandFileName string is equal to the string literal
      * "paramfile", this method returns the same stream as paramFile(). 
      * Otherwise, it returns a reference to a file whose name is given 
      * by the commandFileName string. If setDirectoryId(int) has not 
      * been called, the path to this file (absolute or relative to the 
      * working directory) is equal to the commandFileName string. If 
      * setDirectory() has been called with an integer argument n, the
      * path to this file is obtained adding "n/" as a prefix to the 
      * commandFileName. In either case, if the commandFileName is not
      * "paramfile", the required file is opened for reading the first
      * time it is returned by this method.
      */
      std::istream& commandFile();

      //@}
      
   private:

      /*
      * Name of the parameter file.
      */
      std::string  paramFileName_;

      /*
      * Name of the command file.
      */
      std::string  commandFileName_;

      /*
      * Prefix added to input file names.
      *
      * If this is a file path, with a trailing directory separator,
      * all input files will be sought in this directory.
      */
      std::string  inputPrefix_;

      /*
      * Prefix added to output file names.
      *
      * If this is a file path, with a trailing directory separation, all 
      * output files will be written to this directory.
      */
      std::string  outputPrefix_;

      /*
      * The integer directory id (used in parallel mode).
      */
      std::string  directoryIdPrefix_;

      /*
      * Path for the root directory for this simulation.
      */
      std::string  rootPrefix_;

      /*
      * Pointer to parameter file.
      */
      std::ifstream*  paramFilePtr_;

      /*
      * Pointer to command file.
      */
      std::ifstream*  commandFilePtr_;

      /*
      * Has a directoryId prefix been set?
      */
      bool  hasDirectoryId_;

      /*
      * Has setCommonControl been called?
      */
      bool isCommonControl_;

   };

   /*
   * Return the command file base name.
   */
   inline std::string FileMaster::paramFileName() const
   {  return paramFileName_; }

   /*
   * Return the command file base name.
   */
   inline std::string FileMaster::commandFileName() const
   {  return commandFileName_; }

} 
#endif
