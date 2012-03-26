#ifndef MCMD_MC_DIAGNOSTIC_MANAGER_H
#define MCMD_MC_DIAGNOSTIC_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/DiagnosticManager.h>  // base class

namespace McMd
{

   using namespace Util;

   class Simulation;
   class McSimulation;
   class McSystem;

   /**
   * Manager for Diagnostic objects in an McSimulation.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Diagnostic_Module
   */
   class McDiagnosticManager : public DiagnosticManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      */
      McDiagnosticManager(McSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      * \param system     associated McSystem
      */
      McDiagnosticManager(McSimulation& simulation, McSystem& system);

   protected:

      /**
      * Instantiate and return pointer to a new McDiagnosticFactory object.
      */
      virtual Factory<Diagnostic>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      McSimulation* simulationPtr_;

      /// Pointer to associated McSystem.
      McSystem*   systemPtr_;

   };

}
#endif
