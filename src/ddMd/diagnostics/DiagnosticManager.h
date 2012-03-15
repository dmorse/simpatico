#ifndef DIAGNOSTIC_MANAGER_H
#define DIAGNOSTIC_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnostic.h"                 // template parameter
#include <util/param/Manager.h>         // base class template

namespace DdMd
{

   using namespace Util;

   /**
   * Manager for a list of Diagnostic objects.
   *
   * \ingroup Manager_Module
   * \ingroup Diagnostic_Module
   */
   class DiagnosticManager : public Manager<Diagnostic>
   {

   public:

      /**
      * Constructor.
      */
      DiagnosticManager();

      /**
      * Destructor.
      */
      virtual ~DiagnosticManager();

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      virtual void readParam(std::istream &in);

      /**
      * Call setup method of each Diagnostic.
      * 
      * This method should be called just before the main
      * simulation loop, after an initial configuration is
      * known.
      */
      void setup();
 
      /**
      * Call sample method of each Diagnostic.
      */
      void sample(long iStep);
 
      /**
      * Call output method of each diagnostic.
      */
      void output();
 
   };

}
#endif
