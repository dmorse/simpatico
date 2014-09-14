#ifndef SPAN_CONFIG_READER_H
#define SPAN_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace SpAn
{

   class Configuration;
   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigReader implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup SpAn_ConfigReader_Module
   */
   class ConfigReader  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      ConfigReader();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      */
      ConfigReader(Configuration& configuration, bool needsAuxiliaryFile = false);

      /**
      * Destructor.
      */
      virtual ~ConfigReader();

      /**
      * Read a configuration file.
      *
      * \param file  input file 
      */
      virtual void readConfig(std::ifstream& file) = 0;

      /**
      * Return true if an auxiliary file is needed.
      */
      bool needsAuxiliaryFile() const;

      /**
      * Read auxiliary file.
      * 
      * Empty default implementation.
      *
      * \param file  auxiliary file, open for reading.
      */
      virtual void readAuxiliaryFile(std::ifstream& file)
      {}

   protected:

      Configuration& configuration()
      { return *configurationPtr_; }

   private:

      Configuration* configurationPtr_;

      bool needsAuxiliaryFile_;

   };

   /*
   * Return true if an auxiliary file is needed.
   */
   inline
   bool ConfigReader::needsAuxiliaryFile() const
   {  return needsAuxiliaryFile_; }

}
#endif
