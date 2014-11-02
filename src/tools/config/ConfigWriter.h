#ifndef TOOLS_CONFIG_WRITER_H
#define TOOLS_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace Tools
{

   class Configuration;
   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigWriter implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup Tools_ConfigWriter_Module
   */
   class ConfigWriter  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      ConfigWriter();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      * \param needsAuxiliaryFile Does this class need to read an auxiliary file?
      */
      ConfigWriter(Configuration& configuration, bool needsAuxiliaryFile = false);

      /**
      * Destructor.
      */
      virtual ~ConfigWriter();

      /**
      * Write configuration file.
      *
      * \param file  output file 
      */
      virtual void writeConfig(std::ofstream& file) = 0;

      /**
      * Return true if an auxiliary file is needed.
      */
      bool needsAuxiliaryFile() const;

      /**
      * Read auxiliary file (empty default implementation).
      * 
      * \param file  auxiliary file, open for reading.
      */
      virtual void readAuxiliaryFile(std::ifstream& file)
      {}

   protected:

      Configuration& configuration()
      { return *configurationPtr_; }

   private:

      /// Pointer to parent configuration.
      Configuration* configurationPtr_;

      /// Does this Writer need to read an auxiliary file?
      bool needsAuxiliaryFile_;

   };

   /*
   * Return true if an auxiliary file is needed.
   */
   inline
   bool ConfigWriter::needsAuxiliaryFile() const
   {  return needsAuxiliaryFile_; }

}
#endif
