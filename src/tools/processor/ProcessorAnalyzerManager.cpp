/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ProcessorAnalyzerManager.h" 
#include "ProcessorAnalyzerFactory.h" 

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   ProcessorAnalyzerManager::ProcessorAnalyzerManager(Processor& processor)
   : AnalyzerManager(),
     processorPtr_(&processor)
   {  setClassName("AnalyzerManager"); }

   /*
   * Destructor.
   */
   ProcessorAnalyzerManager::~ProcessorAnalyzerManager()
   {} 

   /*
   * Return pointer to default factory.
   */
   Factory<Analyzer>* ProcessorAnalyzerManager::newDefaultFactory() const
   {
      return new ProcessorAnalyzerFactory(*processorPtr_);
   }
 
}
