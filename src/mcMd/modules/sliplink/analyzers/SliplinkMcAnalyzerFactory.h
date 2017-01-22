#ifndef SLIPLINK_MC_ANALYZER_FACTORY_H
#define SLIPLINK_MC_ANALYZER_FACTORY_H

#include <mcMd/analyzers/mcSystem/McAnalyzerFactory.h>

namespace McMd
{

   /**
   * Custom AnalyzerFactory for an McSimulation
   */
   class SliplinkMcAnalyzerFactory : public McAnalyzerFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SliplinkMcAnalyzerFactory(McSimulation& simulation, McSystem& system)
       : McAnalyzerFactory(simulation, system)
      {}

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

   };

}
#endif
