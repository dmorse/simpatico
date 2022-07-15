/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerMixIn.h"
#include <util/misc/FileMaster.h>

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   AnalyzerMixIn::AnalyzerMixIn(FileMaster& fileMaster) 
    : outputFile_(),
      fileMasterPtr_(&fileMaster)
   {}

   /*
   * Destructor.
   */
   AnalyzerMixIn::~AnalyzerMixIn() 
   {}

   /*
   * Open the output file.
   */
   void AnalyzerMixIn::openOutputFile(std::string name)
   {  fileMasterPtr_->openOutputFile(name, outputFile_); }

}
