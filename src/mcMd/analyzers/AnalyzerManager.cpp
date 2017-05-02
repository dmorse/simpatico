/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerManager.h" 
#include "Analyzer.h" 
#include <util/archives/Serializable_includes.h>

namespace McMd
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
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   void AnalyzerManager::readParameters(std::istream &in)
   {
      read<long>(in,"baseInterval", Analyzer::baseInterval);
      Manager<Analyzer>::readParameters(in);
   }

   /*
   * Call initialize method of each analyzer.
   */
   void AnalyzerManager::setup() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setup();
      }
   }
 
   /*
   * Call sample method of each analyzer.
   */
   void AnalyzerManager::sample(long iStep) 
   {
      UTIL_CHECK(Analyzer::baseInterval > 0);
      UTIL_CHECK(iStep % Analyzer::baseInterval == 0);
      for (int i = 0; i < size(); ++i) {
         (*this)[i].sample(iStep);
      }
   }
 
   /*
   * Call output method of each analyzer.
   */
   void AnalyzerManager::output() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Read instructions for creating objects from file.
   */
   void AnalyzerManager::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<long>(ar, "baseInterval", Analyzer::baseInterval);
      Manager<Analyzer>::loadParameters(ar);
   }

   /*
   * Read instructions for creating objects from file.
   */
   void AnalyzerManager::save(Serializable::OArchive &ar)
   {
      ar << Analyzer::baseInterval;
      Manager<Analyzer>::save(ar);
   }

}
