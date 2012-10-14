#ifndef UTIL_FILE_MASTER_CPP
#define UTIL_FILE_MASTER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "FileMaster.h"
#include <util/global.h>

#include <sstream>

namespace Util
{

   /*
   * Default constructor.
   */
   FileMaster::FileMaster()
    : commandFileName_(),
      inputPrefix_(),
      outputPrefix_(),
      directoryIdPrefix_(),
      rootPrefix_(),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(false),
      isSetParamFileStdIn_(false)
   {  setClassName("FileMaster"); }

   /*
   * Copy constructor.
   */
   FileMaster::FileMaster(const FileMaster& other)
    : inputPrefix_(other.inputPrefix_),
      outputPrefix_(other.outputPrefix_),
      directoryIdPrefix_(other.directoryIdPrefix_),
      rootPrefix_(other.rootPrefix_),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(other.hasDirectoryId_),
      isSetParamFileStdIn_(other.isSetParamFileStdIn_)
   {}

   /*
   * Destructor.
   */
   FileMaster::~FileMaster()
   {
      if (paramFilePtr_) {
         paramFilePtr_->close();
         delete paramFilePtr_;
      }
      if (commandFilePtr_) {
         commandFilePtr_->close();
         delete commandFilePtr_;
      }
   }

   /*
   * Set root prefix for all path names.
   */
   void FileMaster::setRootPrefix(const std::string& rootPrefix)
   {  rootPrefix_ = rootPrefix; }

   /*
   * Set a directory Id prefix.
   */
   void FileMaster::setDirectoryId(int rank)
   {
      std::stringstream sMyId;
      sMyId << rank;
      directoryIdPrefix_ = sMyId.str();
      directoryIdPrefix_ += "/";
      hasDirectoryId_ = true;
   }

   /*
   * Set paramFile() to return std::cin even if a directory id has been set.
   */
   void FileMaster::setParamFileStdIn()
   {  isSetParamFileStdIn_ = true; }

   /*
   * Set the input prefix.
   */
   void FileMaster::setInputPrefix(const std::string& inputPrefix)
   {  inputPrefix_ = inputPrefix; }

   /*
   * Set the output prefix.
   */
   void FileMaster::setOutputPrefix(const std::string& outputPrefix)
   {  outputPrefix_ = outputPrefix; }

   /*
   * Set the output prefix.
   */
   void FileMaster::setCommandFileName(const std::string& commandFileName)
   {  commandFileName_ = commandFileName; }

   /*
   * Read parameters from file.
   */
   void FileMaster::readParameters(std::istream &in)
   {
      read<std::string>(in, "commandFileName",  commandFileName_);
      read<std::string>(in, "inputPrefix",  inputPrefix_);
      read<std::string>(in, "outputPrefix", outputPrefix_);
   }

   #if 0
   /*
   * Load internal state from file.
   */
   void FileMaster::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "commandFileName",  commandFileName_);
      loadParameter<std::string>(ar, "inputPrefix",  inputPrefix_);
      loadParameter<std::string>(ar, "outputPrefix", outputPrefix_);
      ar >> directoryIdPrefix_;
      ar >> rootPrefix_;
      ar >> hasDirectoryId_;
      ar >> isSetParamFileStdIn_;
   }

   /*
   * Save internal state to file.
   */
   void FileMaster::save(Serializable::OArchive &ar)
   {
      ar << commandFileName_;
      ar << inputPrefix_;
      ar << outputPrefix_;
      ar << directoryIdPrefix_;
      ar << rootPrefix_;
      ar << hasDirectoryId_;
      ar << isSetParamFileStdIn_;
   }
   #endif

   /*
   * Get the default parameter stream.
   */
   std::istream& FileMaster::paramFile()
   {
      if ((!hasDirectoryId_) || isSetParamFileStdIn_) {
         return std::cin;
      } else {
         if (!paramFilePtr_) {

            // Construct parameter filename
            std::string filename(rootPrefix_);
            if (hasDirectoryId_) {
               filename += directoryIdPrefix_;
            }
            filename += "param";

            // Open parameter input file
            std::ifstream* filePtr = new std::ifstream();
            filePtr->open(filename.c_str());
            if (filePtr->fail()) {
               std::string message = "Error opening parameter file. Filename: ";
               message += filename;
               UTIL_THROW(message.c_str());
            }

            paramFilePtr_ = filePtr;
         }
         return *paramFilePtr_;
      }
   }

   /*
   * Get the command input stream by reference.
   */
   std::istream& FileMaster::commandFile()
   {
      if (commandFileName_ == "paramfile") {
         return paramFile();
      } else {
         if (!commandFilePtr_) {

            // Construct filename "n/param" for processor n
            std::string filename(rootPrefix_);
            if (hasDirectoryId_) {
               filename += directoryIdPrefix_;
            }
            filename += commandFileName_;

            // Open parameter input file
            std::ifstream* filePtr = new std::ifstream();
            filePtr->open(filename.c_str());
            if (filePtr->fail()) {
               std::string message;
               message = "Error opening command file. Filename: ";
               message += filename;
               UTIL_THROW(message.c_str());
            }
            commandFilePtr_ = filePtr;
         }
         return *commandFilePtr_;
      }
   }

   /*
   * Open and return an input file with specified base name
   */
   void
   FileMaster::openInputFile(const std::string& name, std::ifstream& in) const
   {
      // Construct filename = inputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += inputPrefix_;
      filename += name;

      in.open(filename.c_str());

      // Check for error opening file
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }

   }

   /*
   * Open and return an output file named outputPrefix + name
   */
   void
   FileMaster::openOutputFile(const std::string& name, std::ofstream& out)
   const
   {
      // Construct filename = outputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += outputPrefix_;
      filename += name;

      out.open(filename.c_str());

      // Check for error opening file
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }

   }

   /*
   * Open an input restart parameter file.
   */
   void FileMaster::openParamIFile(const std::string& name, const char* ext,
                                   std::ifstream& in) const
   {
      // Construct filename
      std::string filename(rootPrefix_);
      if (hasDirectoryId_ && !isSetParamFileStdIn_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      filename += ext;

      in.open(filename.c_str());
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open a restart parameter file for writing.
   */
   void FileMaster::openParamOFile(const std::string& name, const char* ext,
                                   std::ofstream& out) const
   {
      // Construct filename 
      std::string filename(rootPrefix_);
      if (hasDirectoryId_ && !isSetParamFileStdIn_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      filename += ext;

      out.open(filename.c_str());
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open an input restart dump file.
   */
   void FileMaster::openRestartIFile(const std::string& name, const char* ext,
                                     std::ifstream& in) const
   {
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      filename += ext;

      in.open(filename.c_str());
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open a restart dump file for writing.
   */
   void FileMaster::openRestartOFile(const std::string& name, const char* ext,
                                     std::ofstream& out) const
   {
      // Construct filename = outputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename = directoryIdPrefix_;
      }
      filename += name;
      filename += ext;

      out.open(filename.c_str());

      // Check for error opening file
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Will paramFile() return std::cin ?
   */
   bool FileMaster::isParamFileStdIn() const
   {  return ((!hasDirectoryId_) || isSetParamFileStdIn_); }

}
#endif
