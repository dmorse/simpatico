#ifndef SLIPLINK_MD_ANALYZER_FACTORY_H
#define SLIPLINK_MD_ANALYZER_FACTORY_H

#include <mcMd/analyzers/mdSystem/MdAnalyzerFactory.h>

namespace McMd
{

   /**
   * Custom AnalyzerFactory for an MdSimulation
   */
   class SliplinkMdAnalyzerFactory : public MdAnalyzerFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SliplinkMdAnalyzerFactory(MdSimulation& simulation, MdSystem& system)
       : MdAnalyzerFactory(simulation, system)
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
