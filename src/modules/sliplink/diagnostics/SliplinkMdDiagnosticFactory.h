#ifndef SLIPLINK_MD_DIAGNOSTIC_FACTORY_H
#define SLIPLINK_MD_DIAGNOSTIC_FACTORY_H

#include <mcMd/diagnostics/mdSystem/MdDiagnosticFactory.h>

namespace McMd
{

   /**
   * Custom DiagnosticFactory for an MdSimulation
   */
   class SliplinkMdDiagnosticFactory : public MdDiagnosticFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SliplinkMdDiagnosticFactory(MdSimulation& simulation, MdSystem& system)
       : MdDiagnosticFactory(simulation, system)
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
