#ifndef TOOLS_CONFIG_WRITER_FACTORY_H
#define TOOLS_CONFIG_WRITER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <tools/config/ConfigWriter.h>

#include <string>

namespace Tools
{

   using namespace Util;
   class Configuration;

   /**
   * Default Factory for subclasses of ConfigWriter.
   *
   * \ingroup Tools_ConfigWriter_Module
   */
   class ConfigWriterFactory : public Factory<ConfigWriter> 
   {

   public:

      /**
      * Constructor
      *
      * \param configuration parent Configuration object
      */
      ConfigWriterFactory(Configuration& configuration);

      /**
      * Create an instance of a specified subclass of ConfigWriter.
      *
      * If the subclassName is recognized, this method returns a
      * pointer to new object. If the name is not recognized, it
      * returns a null pointer.
      *
      * The new object is created using a constructor that takes
      * the parent Tools::Configuration as a parameter. 
      *
      * \param  subclassName  name of a subclass of ConfigWriter 
      * \return ConfigWriter*     pointer to new instance of subclassName
      */
      ConfigWriter* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent Configuration.
      Configuration* configurationPtr_;

   };

}
#endif
