#ifndef SPAN_CONFIG_IO_FACTORY_H
#define SPAN_CONFIG_IO_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <spAn/configIos/ConfigIo.h>

#include <string>

namespace SpAn
{

   using namespace Util;
   class Configuration;

   /**
   * Default Factory for subclasses of ConfigIo.
   *
   * \ingroup SpAn_ConfigIo_Module
   */
   class ConfigIoFactory : public Factory<ConfigIo> 
   {

   public:

      /**
      * Constructor
      *
      * \param configuration parent Configuration object
      */
      ConfigIoFactory(Configuration& configuration);

      /**
      * Create an instance of a specified subclass of ConfigIo.
      *
      * If the subclassName is recognized, this method returns a
      * pointer to new object. If the name is not recognized, it
      * returns a null pointer.
      *
      * The new object is created using a constructor that takes
      * the parent SpAn::Configuration as a parameter. The calling
      * function must invoke DdMd::ConfigIo::initialize() before
      * using the new ConfigIo to read or write configuration
      * files.
      *
      * \param  subclassName  name of a subclass of ConfigIo 
      * \return ConfigIo*     pointer to new instance of subclassName
      */
      ConfigIo* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent Configuration.
      Configuration* configurationPtr_;

   };

}
#endif
