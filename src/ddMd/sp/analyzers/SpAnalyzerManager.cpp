#ifndef DDMD_SP_ANALYZER_MANAGER_CPP
#define DDMD_SP_ANALYZER_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpAnalyzerManager.h" 
#include "SpAnalyzerFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpAnalyzerManager::SpAnalyzerManager(Processor& processor)
   : Manager<SpAnalyzer>(),
     processorPtr_(&processor)
   {  setClassName("SpAnalyzerManager"); }

   /*
   * Destructor.
   */
   SpAnalyzerManager::~SpAnalyzerManager()
   {} 

   /*
   * Read parameter block (without begin and end).
   *
   * \param in input parameter file stream.
   */
   void SpAnalyzerManager::readParameters(std::istream &in)
   {
      Manager<SpAnalyzer>::readParameters(in);
   }

   /*
   * Call setup method of each analyzer.
   */
   void SpAnalyzerManager::setup() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setup();
      }
   }
 
   /*
   * Call clear method of each analyzer.
   */
   void SpAnalyzerManager::clear() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].clear();
      }
   }

   /*
   * Call sample method of each analyzer.
   */
   void SpAnalyzerManager::sample(long iStep) 
   {
      for (int i=0; i < size(); ++i) {
         if ((*this)[i].isAtInterval(iStep)) {
            (*this)[i].sample(iStep);
          }
      }
   }
 
   /*
   * Call output method of each analyzer.
   */
   void SpAnalyzerManager::output() 
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Return pointer to default factory.
   */
   Factory<SpAnalyzer>* SpAnalyzerManager::newDefaultFactory() const
   {
      return new SpAnalyzerFactory(*processorPtr_);
   }
 
}
#endif
