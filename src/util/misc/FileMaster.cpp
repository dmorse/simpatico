#ifndef UTIL_FILE_MASTER_CPP
#define UTIL_FILE_MASTER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
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
    : paramFileName_(),
      commandFileName_(),
      inputPrefix_(),
      outputPrefix_(),
      directoryIdPrefix_(),
      rootPrefix_(),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(false),
      isCommonControl_(false)
   {  setClassName("FileMaster"); }

   /*
   * Copy constructor.
   */
   FileMaster::FileMaster(const FileMaster& other)
    : paramFileName_(),
      commandFileName_(),
      inputPrefix_(other.inputPrefix_),
      outputPrefix_(other.outputPrefix_),
      directoryIdPrefix_(other.directoryIdPrefix_),
      rootPrefix_(other.rootPrefix_),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(other.hasDirectoryId_),
      isCommonControl_(other.isCommonControl_)
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
   void FileMaster::setCommonControl()
   {  isCommonControl_ = true; }

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
      if (commandFileName_.empty()) {
         read<std::string>(in, "commandFileName",  commandFileName_);
      }
      read<std::string>(in, "inputPrefix",  inputPrefix_);
      read<std::string>(in, "outputPrefix", outputPrefix_);
   }

   /*
   * Get the default parameter stream.
   */
   void FileMaster::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "inputPrefix",  inputPrefix_);
      loadParameter<std::string>(ar, "outputPrefix", outputPrefix_);
      ar >> directoryIdPrefix_;
      ar >> rootPrefix_;
      ar >> hasDirectoryId_;
      ar >> isCommonControl_;
      Log::file() << "drectoryIdPrefix_ = " << directoryIdPrefix_ << std::endl;
      Log::file() << "rootPrefix_       = " << rootPrefix_ << std::endl;
      Log::file() << "hasDirectoryId_   = " << hasDirectoryId_ << std::endl;
      Log::file() << "isCommonControl_  = " << isCommonControl_ << std::endl;
   }

   /*
   * Save internal state to file.
   */
   void FileMaster::save(Serializable::OArchive &ar)
   {
      ar << inputPrefix_;
      ar << outputPrefix_;
      ar << directoryIdPrefix_;
      ar << rootPrefix_;
      ar << hasDirectoryId_;
      ar << isCommonControl_;
   }

   /*
   * Get the default parameter stream.
   */
   std::istream& FileMaster::paramFile()
   {
      if (paramFilePtr_) {
         return *paramFilePtr_;
      } else {
         if (paramFileName_.empty()) {
            if (!hasDirectoryId_ || isCommonControl_) {
               return std::cin;
            } else {
               paramFileName_ = "param";
            }
         } 
         paramFilePtr_ = new std::ifstream();
         // Log::file() << "Opening parameter file " 
         //            << paramFileName_  << std::endl;
         openControlFile(paramFileName_, *paramFilePtr_);
         return *paramFilePtr_;
      }
   }

   /*
   * Get the command input stream by reference.
   */
   std::istream& FileMaster::commandFile()
   {
      if (commandFilePtr_) {
         return *commandFilePtr_;
      } else {
         if (commandFileName_.empty()) {
            commandFileName_ = "commands";
         } 
         commandFilePtr_ = new std::ifstream();
         //Log::file() << "Opening command file " 
         //            << paramFileName_  << std::endl;
         openControlFile(commandFileName_, *commandFilePtr_);
         return *commandFilePtr_;
      }
   }

   /*
   * Open an input name with fully specified path.
   */
   void
   FileMaster::open(const std::string& name, std::ifstream& in,
                    std::ios_base::openmode mode) const
   {
      in.open(name.c_str(), mode);
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += name;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open an ofstream with fully specified path.
   */
   void
   FileMaster::open(const std::string& name, std::ofstream& out,
                    std::ios_base::openmode mode) const
   {
      out.open(name.c_str(), mode);
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += name;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open an input restart parameter file.
   */
   void FileMaster::openControlFile(const std::string& name, 
                                   std::ifstream& in) const
   {
      std::string filename(rootPrefix_);
      if (hasDirectoryId_ && !isCommonControl_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      open(filename.c_str(), in);
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
      open(filename.c_str(), in);
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
      open(filename.c_str(), out);
   }

   /*
   * Open and return an input file with specified base name
   */
   void
   FileMaster::openInputFile(const std::string& name, std::ifstream& in,
                             std::ios_base::openmode mode) const
   {
      // Construct filename = inputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += inputPrefix_;
      filename += name;
      open(filename.c_str(), in, mode);
   }

   /*
   * Open and return an input file with specified base name
   */
   void
   FileMaster::openInputFile(const std::string& name, std::ifstream& in) const
   {  openInputFile(name, in, std::ios_base::in); }

   /*
   * Open and return an output file named outputPrefix + name
   */
   void
   FileMaster::openOutputFile(const std::string& name, 
                              std::ofstream& out, 
                              std::ios_base::openmode mode) const
   {
      // Construct filename = outputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += outputPrefix_;
      filename += name;

      open(filename.c_str(), out, mode);
   }

   /*
   * Open and return an output file named outputPrefix + name
   */
   void
   FileMaster::openOutputFile(const std::string& name, 
                              std::ofstream& out, bool append) const
   {
      if (append) {
         openOutputFile(name, out, std::ios::out | std::ios::app);
      } else {
         openOutputFile(name, out, std::ios::out);
      }
   }

   /*
   * Will paramFile() return std::cin ?
   */
   bool FileMaster::isCommonControl() const
   {  return ((!hasDirectoryId_) || isCommonControl_); }

}
#endif
