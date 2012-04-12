#ifndef DDMD_FILE_MASTER_H
#define DDMD_FILE_MASTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>                // base class

#include <fstream>
#include <string>

namespace McMd{ class FileMaster; }

namespace DdMd
{

   using namespace Util;

   /**
   * A FileMaster manages a set of input and output files.
   *
   * A FileMaster manages a set of input and output files in which
   * all input file paths share a common input prefix and all output
   * file paths share a common output prefix.  The openInputFile() 
   * and openOutputFile() methods each open a file using a path that
   * is constructed by prepending the input or output prefix, 
   * respectively, to a filename string parameter. A FileMaster 
   * also provides access to a separate parameter file.
   *
   * The input and output prefix strings are read from file by the
   * readParam() method. To use input files that are all in one 
   * directory and create output files in another, these prefix
   * strings should be set to the desired directory names, followed
   * by a trailing "/" directory separator. For example, to put all
   * input and output files in directories named "in" and "out", the
   * prefix strings should be set to "in/" and "out/", respectively.
   * 
   * The setDirectoryId(int rank) method is designed for use in parallel
   * MPI programs, to allow files associated with different processors 
   * to be put in different directory trees. If this method is called 
   * with a rank parameter n, openInputFile() and openOutputFile() 
   * will thereafter prepend a string "n/" to all file paths before the
   * inputPrefix and output prefix strings. Calling setDirectoryId(n) 
   * and placing prefixes "in/" and "out/" in the parameter file that
   * is read by readParam thus causes the strings "n/in/" and "n/out/" 
   * to be prepended to the names of files opened by openInputFile()
   * and openOutputFile().
   *
   * The paramFile() method returns a default parameter input stream,
   * using different conventions for serial and parallel operation. 
   * If no directory id has been set, as is normally the case for a 
   * serial program, paramFile() returns a reference to the standard
   * input, std::cin.  In parallel programs, if setDirectoryId() has
   * been called, but setParamFileStdIn() has also been called, 
   * paramFile() again returns a reference to std::cin. If a directory
   * id has been set by calling setDirectoryId() has been and 
   * setParamFileStdIn() has not be called, however, paramFile() 
   * returns a reference to a file named "n/param". This file is 
   * opened for reading the first time it is returned by paramFile(), 
   * and is closed by the FileMaster destructor.
   *
   */
   class FileMaster : public ParamComposite
   {
   
   public:

      /**
      * Constructor.
      */
      FileMaster();
 
      /**
      * Constructor.
      */
      FileMaster(const McMd::FileMaster& other);
 
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
      virtual void readParam(std::istream& in);

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
      * Open an input file named [directoryIdPrefix] + inputPrefix + filename.
      *
      * The directoryIdPrefix "n/" is added as a prefix to the path name if and 
      * only if setDirectoryId(int) has been called with an argument n. 
      *
      * \param  filename file name, without any prefix
      * \param  in       ifstream object to associated with a file
      */
      void openInputFile(const std::string& filename, std::ifstream& in) const;

      /**
      * Open an output file named [directoryIdPrefix] + inputPrefix + filename.
      *
      * The directoryIdPrefix "n/" is added as a prefix to the path name if and 
      * only if setDirectoryId(int) has been called with an argument n. 
      *
      * \param  filename  file name, without any prefix
      * \param  out       ofstream object to associated with a file
      */
      void 
      openOutputFile(const std::string& filename, std::ofstream& out) const;

      /**
      * Open an input restart parameter file.
      *
      * \param  name  base file name, without any prefix
      * \param  ext   file name extension.
      * \param  in    ifstream object to open
      */
      void openParamIFile(const std::string& name, const char* ext,
                          std::ifstream& in) const;

      /**
      * Open a restart parameter file for writing.
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  out   ofstream object to open
      */
      void openParamOFile(const std::string& name, const char* ext,
                          std::ofstream& out) const;

      /**
      * Open an input restart dump file.
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  in    ifstream object to open
      */
      void openRestartIFile(const std::string& name, const char* ext,
                            std::ifstream& in) const;

      /**
      * Open a restart dump file for writing.
      *
      * \param  name  base file name, without any prefix or extension
      * \param  ext   file name extension.
      * \param  out   ofstream object to open
      */
      void openRestartOFile(const std::string& name, const char* ext, 
                            std::ofstream& out) const;

      /**
      * Return the command file name.
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
