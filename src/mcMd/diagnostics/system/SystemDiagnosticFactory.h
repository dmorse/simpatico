#ifndef SYSTEM_DIAGNOSTIC_FACTORY_H
#define SYSTEM_DIAGNOSTIC_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <mcMd/diagnostics/Diagnostic.h>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;

   /**
   * DiagnosticFactory for any System (for mc or md).
   *
   * \ingroup Factory_Module
   * \ingroup Diagnostic_Module
   */
   class SystemDiagnosticFactory : public Factory<Diagnostic>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SystemDiagnosticFactory(Simulation& simulation, System& system);

      /** 
      * Return pointer to a new Diagnostic object.
      *
      * \param  className name of a subclass of Diagnostic.
      * \return base class pointer to a new instance of className.
      */
      virtual Diagnostic* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Simulation.
      */
      Simulation& simulation() const
      {  return *simulationPtr_; }

   private:

      /**
      * Return reference to parent System.
      */
      System& system() const
      {  return *systemPtr_; }

      Simulation *simulationPtr_;
      System     *systemPtr_;

   };

}
#endif
