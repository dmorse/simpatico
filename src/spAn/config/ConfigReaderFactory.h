#ifndef SPAN_CONFIG_READER_FACTORY_H
#define SPAN_CONFIG_READER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <spAn/config/ConfigReader.h>

#include <string>

namespace SpAn
{

   using namespace Util;
   class Configuration;

   /**
   * Default Factory for subclasses of ConfigReader.
   *
   * \ingroup SpAn_ConfigReader_Module
   */
   class ConfigReaderFactory : public Factory<ConfigReader> 
   {

   public:

      /**
      * Constructor
      *
      * \param configuration parent Configuration object
      */
      ConfigReaderFactory(Configuration& configuration);

      /**
      * Create an instance of a specified subclass of ConfigReader.
      *
      * If the subclassName is recognized, this method returns a
      * pointer to new object. If the name is not recognized, it
      * returns a null pointer.
      *
      * The new object is created using a constructor that takes
      * the parent SpAn::Configuration as a parameter. The calling
      * function must invoke DdMd::ConfigReader::initialize() before
      * using the new ConfigReader to read or write configuration
      * files.
      *
      * \param  subclassName  name of a subclass of ConfigReader 
      * \return ConfigReader*     pointer to new instance of subclassName
      */
      ConfigReader* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent Configuration.
      Configuration* configurationPtr_;

   };

}
#endif
