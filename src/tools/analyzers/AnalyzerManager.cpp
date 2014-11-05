/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerManager.h" 

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   AnalyzerManager::AnalyzerManager()
   : Manager<Analyzer>()
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
      Manager<Analyzer>::readParameters(in);
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
      for (int i=0; i < size(); ++i) {
         if ((*this)[i].isAtInterval(iStep)) {
            (*this)[i].sample(iStep);
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

}
