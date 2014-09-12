#ifndef SPAN_CONFIG_WRITER_CPP
#define SPAN_CONFIG_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigWriter.h"

namespace SpAn
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
   ConfigWriter::ConfigWriter(Configuration& configuration)
    : configurationPtr_(&configuration)
   {  setClassName("ConfigWriter"); }

   /*
   * Destructor.
   */
   ConfigWriter::~ConfigWriter()
   {}

}
#endif
