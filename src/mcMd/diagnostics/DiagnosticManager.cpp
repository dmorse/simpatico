#ifndef MCMD_DIAGNOSTIC_MANAGER_CPP
#define MCMD_DIAGNOSTIC_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DiagnosticManager.h" 
#include "Diagnostic.h" 
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DiagnosticManager::DiagnosticManager()
   : Manager<Diagnostic>()
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
   * Call initialize method of each diagnostic.
   */
   void DiagnosticManager::initialize() 
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].initialize();
      }
   }
 
   /*
   * Call sample method of each diagnostic.
   */
   void DiagnosticManager::sample(long iStep) 
   {
      if (iStep % Diagnostic::baseInterval == 0) { 
         for (int i=0; i < size(); ++i) {
            (*this)[i].sample(iStep);
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
   * Save state to binary file archive.
   */
   void DiagnosticManager::save(Serializable::OArchiveType& ar)
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].save(ar);
      }
   }

   /*
   * Load state from a binary file archive.
   */
   void DiagnosticManager::load(Serializable::IArchiveType& ar)
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].load(ar);
      }
   }

}
#endif
