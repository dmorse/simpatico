#ifndef MCMD_SYSTEM_ANALYZER_FACTORY_H
#define MCMD_SYSTEM_ANALYZER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <mcMd/analyzers/Analyzer.h>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;

   /**
   * AnalyzerFactory for any System (for mc or md).
   *
   * \ingroup McMd_Factory_Module
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class SystemAnalyzerFactory : public Factory<Analyzer>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SystemAnalyzerFactory(Simulation& simulation, System& system);

      /** 
      * Return pointer to a new Analyzer object.
      *
      * \param  className name of a subclass of Analyzer.
      * \return base class pointer to a new instance of className.
      */
      virtual Analyzer* factory(const std::string& className) const;

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
