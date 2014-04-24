#ifndef MDPP_CONFIG_IO_H
#define MDPP_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace MdPp
{

   class Processor;
   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigIo implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup MdPp_ConfigIo_Module
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
      * \param simulation parent Simulation object.
      */
      ConfigIo(Processor& processor);

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

   private:

      Processor* processorPtr_;

   };

}
#endif
