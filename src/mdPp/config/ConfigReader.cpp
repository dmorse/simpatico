/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigReader.h"
#include <mdPp/storage/Configuration.h>

#include <iostream>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigReader::ConfigReader()
     : configurationPtr_(0)
   {  setClassName("ConfigReader"); }

   /*
   * Constructor.
   */
   ConfigReader::ConfigReader(Configuration& configuration, 
                              bool needsAuxiliaryFile)
     : configurationPtr_(&configuration),
       needsAuxiliaryFile_(needsAuxiliaryFile)
   {  setClassName("ConfigReader"); }

   /*
   * Destructor.
   */
   ConfigReader::~ConfigReader()
   {}

}
