#ifndef MCMD_MD_DIAGNOSTIC_MANAGER_H
#define MCMD_MD_DIAGNOSTIC_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/DiagnosticManager.h>  // base class

namespace Util { template <typename T> class Factory; }

namespace McMd
{

   using namespace Util;

   class Simulation;
   class MdSimulation;
   class MdSystem;

   /**
   * Manager for Diagnostic objects in an MdSimulation.
   *
   * The default Factory<Diagnostic> object for an MdDiagnosticManager is
   * an MdDiagnosticFactory.
   *
   * \ingroup Manager_Module
   * \ingroup Diagnostic_Module
   */
   class MdDiagnosticManager : public DiagnosticManager
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent MdSimulation
      */
      MdDiagnosticManager(MdSimulation& simulation);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation
      * \param system     associated MdSystem
      */
      MdDiagnosticManager(MdSimulation& simulation, MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~MdDiagnosticManager();

   protected:

      /**
      * Create and return pointer to a new MdDiagnosticFactory object
      */
      virtual Util::Factory<Diagnostic>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation
      MdSimulation* simulationPtr_;

      /// Pointer to associated MdSystem.
      MdSystem*   systemPtr_;

   };

}
#endif
