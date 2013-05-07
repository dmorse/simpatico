#ifndef MCMD_MC_DIAGNOSTIC_FACTORY_H
#define MCMD_MC_DIAGNOSTIC_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>     // base class template
#include <mcMd/diagnostics/Diagnostic.h>    // base class template parameter
#include <mcMd/diagnostics/system/SystemDiagnosticFactory.h>    // member

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * DiagnosticFactory for an McSimulation
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Diagnostic_Module
   */
   class McDiagnosticFactory : public Factory<Diagnostic>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      * \param system     parent McSystem
      */
      McDiagnosticFactory(McSimulation& simulation, McSystem& system);

      /** 
      * Return pointer to a new Diagnostic object.
      *
      * \param  className name of a subclass of Diagnostic.
      * \return base class pointer to a new instance of className.
      */
      virtual Diagnostic* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent McSystem.
      */
      McSystem& system() const
      { return *systemPtr_; }

      /**
      * Return reference to parent McSimulation.
      */
      McSimulation& simulation() const
      { return *simulationPtr_; }

   private:

      // Factory for diagnostics for any System
      SystemDiagnosticFactory systemFactory_;

      // Pointer to parent Simulation
      McSimulation *simulationPtr_;

      // Pointer to parent System
      McSystem  *systemPtr_;

   };

}
#endif
