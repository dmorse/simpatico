#ifndef UTIL_FILE_MASTER_H
#define UTIL_FILE_MASTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>      // base class

#include <fstream>
#include <string>

namespace Util
{

   /**
   * A FileMaster manages a set of input and output files.
   *
   * A FileMaster manages a set of input and output files in which
   * all input file paths share a common input prefix and all output
   * file paths share a common output prefix.  The openInputFile() 
   * and openOutputFile() methods each open a file using a path 
   * that is constructed by a prefix to a base filename that is
   * passed to the function as a parameter, using different prefixes
   * for input and output files.  A FileMaster also provides access 
   * to a separate parameter file.
   *
   * The paths constructed by a FileMaster are made by concatenating
   * elements that may include any of the following (presented in the
   * order they appear:
   *
   *  - An optional root prefix, which is the path from the current
   *    working directory to the root level directory for input and 
   *    output files for this run. The root prefix defaults to an 
   *    empty string, and may be set by the setRootPrefix method. 
   *    If it is set, the root prefix should end with a "/" directory 
   *    separator. 
   *
   *  - In MPI simulations of multiple systems, filenames may contain
   *    a directory id prefix of the form "p/", where p is an integer
   *    id that identifies a system. The integer directory id may be
   *    set with the setDirectoryId(int id) method. A directory id 
   *    is added to paths only if a directory id has been set. 
   *
   *  - Input files opened by the openInputFile() method and output 
   *    files opened by the openOutputFile() method contain an 
   *    inputPrefix or outputPrefix string, respectively. The 
   *    inputPrefix and outputPrefix strings are read from the
   *    parameter file by the readParam method.
   *
   * The paths to input and output files opened by the openInputFile()
   * and openOutputFile() methods are constructed by concatenating the
   * optional rootPrefix (if any), the optional directory id string
   * (if any), the inputPrefix or outputPrefix string, and a base name
   * that is passed explicitly to the method. 
   *
   * To put all input and output files for a given system in different
   * directories, the input and output prefix strings should be set to
   * the names of the desired directory names, followed by a trailing
   * "/" directory separator.  For example, to put all input and 
   * output files in directories named "in" and "out", these prefix 
   * strings should be set to "in/" and "out/", respectively.
   * 
   * The setDirectoryId(int rank) method is designed for use in parallel
   * MPI programs, to allow files associated with different systems to
   * be put in different directory trees. If this method is called with
   * an 3, while using input and output prefixes "in/" and "out/" and 
   * no overall rootPrefix, the openInputFile() and openOutputFile() 
   * functions will open all input and output files in directories 
   * "3/in" and "3/out", respectively.
   *
   * The paramFile() method returns a default parameter input stream.
   * If no directory id has been set, as is normally the case for
   * simulations of a single system, paramFile() returns a reference 
   * to the standard input, std::cin.  If a directory id has been set
   * and setParamFileStdIn() has also been called, paramFile() again 
   * returns a reference to std::cin. If a directory id has been set 
   * and setParamFileStdIn() has not be called, however, paramFile() 
   * returns a reference to a file named "n/param" for a processor
   * with directoryId == n. If a rootPrefix has been set, the root 
   * prefix is prepended to the path "n/param". This file is opened
   * for reading the first time it is returned by paramFile(), and 
   * is closed by the FileMaster destructor.
   *
   * The class also provides methods to open a command file, and 
   * parameter and restart input and output files that are used for
   * restarting a simulation. The paths to each of these files are
   * constructed as a concatenation of a root prefix (if any), a
   * directory id string (if any), and a base name provided as a
   * parameter to these methods.
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
      * Open an output file.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + outputPrefix + filename.
      *
      * \param  filename  file name, without any prefix
      * \param  out       ofstream object to associated with a file
      * \param  append    True if we are appending to the file (default=false)
      */
      void 
      openOutputFile(const std::string& filename, std::ofstream& out, bool append = false) const;

      /**
      * Open an input restart parameter file for reading.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id is set
      * and the setParamFileStdIn method has not been set. 
      *
      * \param  name  base file name, without any prefix
      * \param  ext   file name extension.
      * \param  in    ifstream object to open
      */
      void openParamIFile(const std::string& name, const char* ext,
                          std::ifstream& in) const;

      /**
      * Open an output restart parameter file for writing.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * The directoryIdPrefix is included only if a directory id is set
      * and the setParamFileStdIn method has not been set. 
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  out   ofstream object to open
      */
      void openParamOFile(const std::string& name, const char* ext,
                          std::ofstream& out) const;

      /**
      * Open an input restart dump file for reading.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  in    ifstream object to open
      */
      void openRestartIFile(const std::string& name, const char* ext,
                            std::ifstream& in) const;

      /**
      * Open an output restart dump file for writing.
      *
      * The path to this file constructed by concatenating:
      * [rootPrefix] + [directoryIdPrefix] + name + "." + ext
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  out   ofstream object to open
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
      std::string  commandFileName_;

      /*
      * Prefix added to input file names.
      *
      * If this is a file path, with a trailing directory separator,
      * all input files will be sought in this directory.
      */
      std::string    inputPrefix_;

      /*
      * Prefix added to output file names.
      *
      * If this is a file path, with a trailing directory separation, all 
      * output files will be written to this directory.
      */
      std::string    outputPrefix_;

      /*
      * The integer directory id (used in parallel mode).
      */
      std::string    directoryIdPrefix_;

      /*
      * Path for the root directory for this simulation.
      */
      std::string    rootPrefix_;

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
