#ifndef DDMD_DIAGNOSTIC_MANAGER_H
#define DDMD_DIAGNOSTIC_MANAGER_H

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

   class Simulation;

   /**
   * Manager for a list of Diagnostic objects.
   *
   * \ingroup DdMd_Manager_Module
   * \ingroup DdMd_Diagnostic_Module
   */
   class DiagnosticManager : public Manager<Diagnostic>
   {

   public:

      /**
      * Constructor.
      */
      DiagnosticManager(Simulation& simulation);

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

      /**
      * Return pointer to a new default factory.
      */
      virtual Factory<Diagnostic>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;
 
   };

}
#endif
