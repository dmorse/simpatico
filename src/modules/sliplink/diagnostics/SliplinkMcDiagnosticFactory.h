#ifndef SLIPLINK_MC_DIAGNOSTIC_FACTORY_H
#define SLIPLINK_MC_DIAGNOSTIC_FACTORY_H

#include <mcMd/diagnostics/mcSystem/McDiagnosticFactory.h>

namespace McMd
{

   /**
   * Custom DiagnosticFactory for an McSimulation
   */
   class SliplinkMcDiagnosticFactory : public McDiagnosticFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SliplinkMcDiagnosticFactory(McSimulation& simulation, McSystem& system)
       : McDiagnosticFactory(simulation, system)
      {}

      /** 
      * Return pointer to a new Diagnostic object.
      *
      * \param  className name of a subclass of Diagnostic.
      * \return base class pointer to a new instance of className.
      */
      virtual Diagnostic* factory(const std::string& className) const;

   };

}
#endif
