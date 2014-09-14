#ifndef SPAN_CONFIG_WRITER_FACTORY_H
#define SPAN_CONFIG_WRITER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <spAn/config/ConfigWriter.h>

#include <string>

namespace SpAn
{

   using namespace Util;
   class Configuration;

   /**
   * Default Factory for subclasses of ConfigWriter.
   *
   * \ingroup SpAn_ConfigWriter_Module
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
      * the parent SpAn::Configuration as a parameter. 
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
