#ifndef MDPP_CONFIG_READER_H
#define MDPP_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace Simp {
   template <int N> class SpeciesGroup;
}

namespace MdPp
{

   class Configuration;
   template <int N> class GroupStorage;
   using namespace Util;
   using namespace Simp;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigReader implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup MdPp_ConfigReader_Module
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
      * \param needsAuxiliaryFile Does this class need an auxiliary input file?
      */
      ConfigReader(Configuration& configuration, 
                   bool needsAuxiliaryFile = false);

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

      /**
      * Get parent Configuration by reference.
      */
      Configuration& configuration()
      { return *configurationPtr_; }

   private:

      /// Pointer to parent Configuration.
      Configuration* configurationPtr_;

      /// Does this reader require an auxiliary file?
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
