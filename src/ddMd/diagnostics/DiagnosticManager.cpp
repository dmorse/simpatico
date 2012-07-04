#ifndef DDMD_DIAGNOSTIC_MANAGER_CPP
#define DDMD_DIAGNOSTIC_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DiagnosticManager.h" 
#include "DiagnosticFactory.h" 

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DiagnosticManager::DiagnosticManager(Simulation& simulation)
   : Manager<Diagnostic>(),
     simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   DiagnosticManager::~DiagnosticManager()
   {} 

   /*
   * Read parameter file. 
   *
   * \param in input parameter file stream.
   */
   void DiagnosticManager::readParam(std::istream &in)
   {
      readBegin(in, "DiagnosticManager");
      read<long>(in,"baseInterval", Diagnostic::baseInterval);
      Manager<Diagnostic>::readParam(in);
   }

   /*
   * Call setup method of each diagnostic.
   */
   void DiagnosticManager::setup() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].setup();
      }
   }
 
   /*
   * Call sample method of each diagnostic.
   */
   void DiagnosticManager::sample(long iStep) 
   {
      if (iStep % Diagnostic::baseInterval == 0) { 
         for (int i=0; i < size(); ++i) {
            if ((*this)[i].isAtInterval(iStep)) {
               (*this)[i].sample(iStep);
            }
         }
      }
   }
 
   /*
   * Call output method of each diagnostic.
   */
   void DiagnosticManager::output() 
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Return pointer to default factory.
   */
   Factory<Diagnostic>* DiagnosticManager::newDefaultFactory() const
   {
      return new DiagnosticFactory(*simulationPtr_);
   }
 
}
#endif
