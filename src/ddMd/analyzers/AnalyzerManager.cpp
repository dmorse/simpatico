#ifndef DDMD_ANALYZER_MANAGER_CPP
#define DDMD_ANALYZER_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerManager.h" 
#include "AnalyzerFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AnalyzerManager::AnalyzerManager(Simulation& simulation)
   : Manager<Analyzer>(),
     simulationPtr_(&simulation)
   {  setClassName("AnalyzerManager"); }

   /*
   * Destructor.
   */
   AnalyzerManager::~AnalyzerManager()
   {} 

   /*
   * Read parameter block (without begin and end).
   *
   * \param in input parameter file stream.
   */
   void AnalyzerManager::readParameters(std::istream &in)
   {
      read<long>(in,"baseInterval", Analyzer::baseInterval);
      Manager<Analyzer>::readParameters(in);
   }

   /*
   * Load internal state from an archive.
   */
   void AnalyzerManager::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<long>(ar, "baseInterval", Analyzer::baseInterval);
      Manager<Analyzer>::loadParameters(ar);
   }

   /*
   * Save internal state to an archive.
   */
   void AnalyzerManager::save(Serializable::OArchive &ar)
   {
      ar << Analyzer::baseInterval;
      Manager<Analyzer>::save(ar);
   }

  
   /*
   * Call setup method of each analyzer.
   */
   void AnalyzerManager::setup() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setup();
      }
   }
 
   /*
   * Call clear method of each analyzer.
   */
   void AnalyzerManager::clear() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].clear();
      }
   }

   /*
   * Call sample method of each analyzer.
   */
   void AnalyzerManager::sample(long iStep) 
   {
      if (Analyzer::baseInterval > 0) {
         if (iStep % Analyzer::baseInterval == 0) { 
            for (int i=0; i < size(); ++i) {
               if ((*this)[i].isAtInterval(iStep)) {
                  (*this)[i].sample(iStep);
               }
            }
         }
      }
   }
 
   /*
   * Call output method of each analyzer.
   */
   void AnalyzerManager::output() 
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Return pointer to default factory.
   */
   Factory<Analyzer>* AnalyzerManager::newDefaultFactory() const
   {
      return new AnalyzerFactory(*simulationPtr_);
   }
 
}
#endif
