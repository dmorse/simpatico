#ifndef DDMD_SP_CONFIG_IO_H
#define DDMD_SP_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace DdMd
{

   class SpConfiguration;
   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of SpConfigIo implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup DdMd_SpConfigIo_Module
   */
   class SpConfigIo  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      SpConfigIo();

      /**
      * Constructor.
      *
      * \param configuration parent SpConfiguration object
      */
      SpConfigIo(SpConfiguration& configuration);

      /**
      * Destructor.
      */
      virtual ~SpConfigIo();

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

      SpConfiguration& configuration()
      { return *configurationPtr_; }

   private:

      SpConfiguration* configurationPtr_;

   };

}
#endif
