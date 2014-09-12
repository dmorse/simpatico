#ifndef SPAN_CONFIG_WRITER_H
#define SPAN_CONFIG_WRITER_H

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
   * Each concrete subclass of ConfigWriter implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup SpAn_ConfigWriter_Module
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
      */
      ConfigWriter(Configuration& configuration);

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

   protected:

      Configuration& configuration()
      { return *configurationPtr_; }

   private:

      Configuration* configurationPtr_;

   };

}
#endif
