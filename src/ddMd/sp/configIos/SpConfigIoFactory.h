#ifndef DDMD_SP_CONFIG_IO_FACTORY_H
#define DDMD_SP_CONFIG_IO_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <ddMd/sp/configIos/SpConfigIo.h>

#include <string>

namespace DdMd
{

   using namespace Util;
   class SpConfiguration;

   /**
   * Default Factory for subclasses of SpConfigIo.
   *
   * \ingroup DdMd_SpConfigIo_Module
   */
   class SpConfigIoFactory : public Factory<SpConfigIo> 
   {

   public:

      /**
      * Constructor
      *
      * \param configuration parent SpConfiguration object
      */
      SpConfigIoFactory(SpConfiguration& configuration);

      /**
      * Create an instance of a specified subclass of SpConfigIo.
      *
      * If the subclassName is recognized, this method returns a
      * pointer to new object. If the name is not recognized, it
      * returns a null pointer.
      *
      * The new object is created using a constructor that takes
      * the parent DdMd::Simulation as a parameter. The calling
      * function must invoke DdMd::SpConfigIo::initialize() before
      * using the new SpConfigIo to read or write configuration
      * files.
      *
      * \param  subclassName  name of a subclass of SpConfigIo 
      * \return SpConfigIo*     pointer to new instance of subclassName
      */
      SpConfigIo* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent SpConfiguration.
      SpConfiguration* configurationPtr_;

   };

}
#endif
