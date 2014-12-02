/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigWriter.h"

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigWriter::ConfigWriter()
     : configurationPtr_(0)
   {  setClassName("ConfigWriter"); }

   /*
   * Constructor.
   */
   ConfigWriter::ConfigWriter(Configuration& configuration, bool needsAuxiliaryFile)
    : configurationPtr_(&configuration),
      needsAuxiliaryFile_(needsAuxiliaryFile)
   {  setClassName("ConfigWriter"); }

   /*
   * Destructor.
   */
   ConfigWriter::~ConfigWriter()
   {}

}
