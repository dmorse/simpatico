#ifndef SPAN_CONFIG_IO_H
#define SPAN_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Each concrete subclass of ConfigIo implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup SpAn_ConfigIo_Module
   */
   class ConfigIo  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      ConfigIo();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      */
      ConfigIo(Configuration& configuration);

      /**
      * Destructor.
      */
      virtual ~ConfigIo();

      /**
      * Read a configuration file.
      *
      * \param file  input file 
      */
      virtual void readConfig(std::ifstream& file) = 0;

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
