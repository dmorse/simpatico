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
   * A FileMaster manages a set of input and output files.
   *
   * A FileMaster manages provides methods to open input and output 
   * files within several locations within a standard directory 
   * structure for molecular simulations. The directory structure 
   * is slightly different for simulations of a single physical 
   * system then for simulations of multiple independent systems, 
   * as discussed below.
   *  
   * For simulations of a single physical system (whether 
   * performed by either a serial or parallel program), we define 
   * three prefix strings for path names, which are usually 
   * correspond to different directories within a tree:
   *
   *   - an root directory path prefix
   *   - a input data path prefix
   *   - an output data path prefix
   *
   * Paths to command and restart files are constructed by
   * pre-pending the root directory path prefix to base names to 
   * base name for these files.  The root directory path prefix must
   * be the path to the root directory for a simulation from the current 
   * working directory. It should either be an empty string (to indicate
   * the current working directory) or a valid directory name that
   * ends with a directory separator "/". The root prefix is set to 
   * an empty string by default. Paths to input and output data files
   * files are constructed by prepending first the root directory 
   * path prefix, then the input or output prefix to base names for 
   * specific data files. The input and output prefixes are often
   * chosen to be names of subdirectories of the root directory 
   * that hold input and output data files (such as "in/" and "out/"), 
   * in which case they should end with a directory separator.
   *
   * For parallel simulations of a multiple independent systems, the
   * root directory must contain a set of directories whose names are 
   * integers corresponding to integer indices for different systems. 
   * Thus for example, the directory that contains files that are 
   * specific to system number 6 would be a subdirectory "6/" of the 
   * root directory. In this case, path names for different types of 
   * files are thus constructed from four path prefixes:
   *
   *   - the root directory path prefix
   *   - an integer directory id prefix, such as "6/" 
   *   - an input data path prefix 
   *   - an output data path prefix 
   *
   * In this case, the path to an input or output file containing 
   * data for a particular system is constructed by concatenting
   * the root directory path prefix, the directory id prefix, the
   * input or output prefix, and a base name for a specific file.
   * When the input and output prefixes are directory paths,
   * such as "in/" and "out/", data files associated with different 
   * systems in separate subdirectories of the associated numbered 
   * system level directories. In multi-system simulations, restart 
   * files for a specific system are placed in the associated 
   * numbered system directory.
   *
   * The location of parameter and command files in multi-system
   * simulations depends on whether the setParamStdIn() method 
   * has been called. If this method has not been called, separate
   * parameter and command files for each system are read from 
   * the numbered system directory for that system. If this method
   * has been called, a single parameter file is read from standard
   * in and a single command file is read from the root directory,
   * using a path that is constructed by concatenating the root
   * prefix with a base name for the command file. 
   *
   * The root prefix and system id prefix are empty strings by
   * default, but may be set by calling the setRootPrefix() and
   * setDirectoryId() member functions, respectively. 
   * 
   * The openInputFile() and openOutputFile() functions open input
   * and output data files, using the above conventions for paths.
   * The openRestartIFile() and openRestartOFile() open input and
   * output restart files, respectively. 
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

      /**
      * Set the path from current directory to root data directory.
      *
      * \param rootPrefix root prefix for all file names.
      */
      void setRootPrefix(const std::string& rootPrefix);

      /**
      * Set an integer directory identifier for this processor.
      *
      * This method should be called only for parallel operation. The
      * directoryId is normally set to be the MPI rank of this processor.
      * After calling setDirectoryId(n) with an integer n, a prefix "n/" 
      * will be prepended to the paths of all input and output file 
      * files that are opened by openInputFile() and openOutputFile().
      *
      * \param directoryId integer directory name for input and output files.
      */
      void setDirectoryId(int directoryId);

      /**
      * Set paramFile() to return std::cin even if a directory id is set.
      */ 
      void setParamFileStdIn();

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
      * Set the command file name.
      *
      * \param commandFileName name of command file
      */
      void setCommandFileName(const std::string& commandFileName);
 
      /**
      * Read parameter file.
      *
      * Reads the inputPrefix and outputPrefix string variables.
      *
      * \param in pararameter file input stream
      */
      virtual void readParameters(std::istream& in);

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
   
      /**
      * Get a default parameter stream by reference.
      *
      * If setDirectoryId() has not been called, of if setParamFileStdIn()
      * has been called, this method returns std::cin.
      * 
      * If setDirectoryId() has been called and setParamFileStdInd() has not,
      * this method returns a reference to a file "n/param". This file is 
      * opened for reading the first time it is returned by this method.
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

      /**
      * Open an input file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + inputPrefix + filename.
      *
      * \param  filename file name, without any prefix
      * \param  in       ifstream object to associated with a file
      */
      void openInputFile(const std::string& filename, std::ifstream& in) const;

      /**
      * Open an input file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + inputPrefix + filename.
      *
      * \param  filename file name, without any prefix
      * \param  in       ifstream object to associated with a file
      * \param  mode     bit mask that specifies opening mode
      */
      void openInputFile(const std::string& filename, std::ifstream& in, 
                         std::ios_base::openmode mode) const;

      /**
      * Open an output file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + outputPrefix + filename.
      *
      * \param  filename  file name, without any prefix
      * \param  out       ofstream object to associated with a file
      * \param  append    True if appending to the file (default=false)
      */
      void 
      openOutputFile(const std::string& filename, std::ofstream& out, 
                     bool append = false) const;

      /**
      * Open an output file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + outputPrefix + filename.
      *
      * The directoryIdPrefix is included only if a directory id 
      * is set and setParamFileStdIn() has not been invoked.
      *
      * \param  filename  file name, without any prefix
      * \param  out       ofstream object to associated with a file
      * \param  mode      bit mask that specifies opening mode
      */
      void 
      openOutputFile(const std::string& filename, std::ofstream& out, 
                     std::ios_base::openmode mode) const;
 
      /**
      * Open an input file in the parameter directory.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id 
      * is set and setParamFileStdIn() has not been invoked.
      *
      *
      * \param  name  base file name, without any prefix
      * \param  ext   file name extension.
      * \param  in    ifstream object to open
      */
      void openParamIFile(const std::string& name, const char* ext,
                          std::ifstream& in) const;

      /**
      * Open an input file in the parameter directory.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id 
      * is set and setParamFileStdIn() has not been invoked.
      *
      *
      * \param name  base file name, without any prefix or extension
      * \param ext  file name extension.
      * \param out  ofstream object to open
      */
      void openParamOFile(const std::string& name, const char* ext,
                          std::ofstream& out) const;

      /**
      * Open an input restart file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id 
      * is set and the setParamFileStdIn method has not been set. 
      *
      * \param name  base file name, without any prefix or extension
      * \param ext  file name extension.
      * \param in  ifstream object to open
      */
      void openRestartIFile(const std::string& name, const char* ext,
                            std::ifstream& in) const;

      /**
      * Open an output restart file. 
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * \param name  base file name, without any prefix or extension
      * \param ext  file name extension.
      * \param out  ofstream object to open
      */
      void openRestartOFile(const std::string& name, const char* ext, 
                            std::ofstream& out) const;

      /**
      * Return the command file name.
      *
      * The base name of the param file is read from the parameter file
      * by the readParameters() method. 
      */
      std::string commandFileName() const;

      /**
      * Is paramFile() set to return std::cin ?
      */
      bool isParamFileStdIn() const;

   private:

      /*
      * Name of the command file.
      */
      std::string commandFileName_;

      /*
      * Prefix added to input file names.
      *
      * If this is a file path, with a trailing directory separator,
      * all input files will be sought in this directory.
      */
      std::string inputPrefix_;

      /*
      * Prefix added to output file names.
      *
      * If this is a file path, with a trailing directory separation, all 
      * output files will be written to this directory.
      */
      std::string outputPrefix_;

      /*
      * The integer directory id (used in parallel mode).
      */
      std::string directoryIdPrefix_;

      /*
      * Path for the root directory for this simulation.
      */
      std::string rootPrefix_;

      /*
      * Pointer to parameter file.
      */
      std::ifstream* paramFilePtr_;

      /*
      * Pointer to command file.
      */
      std::ifstream* commandFilePtr_;

      /*
      * Has a directoryId prefix been set?
      */
      bool hasDirectoryId_;

      /*
      * Has setParamFileStdIn been called?
      */
      bool isSetParamFileStdIn_;

   };

   /*
   * Return the command file base name.
   */
   inline std::string FileMaster::commandFileName() const
   { return commandFileName_; }

} 
#endif
