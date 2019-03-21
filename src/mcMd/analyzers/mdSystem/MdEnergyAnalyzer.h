#ifndef MCMD_MD_ENERGY_ANALYZER_H
#define MCMD_MD_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageListAnalyzer.h>
#include <mcMd/mdSimulation/MdSystem.h>

namespace Util{
   class Average;
}

namespace McMd
{

   using namespace Util;

   /**
   * Compute averages and output block averages of energy components.
   *
   * This class computes separate averages for each component of the
   * total simulation energy (kinetic, pair, bond, etc.) as well as
   * for the total, and periodically outputs block averages of each
   * to a file.
   *
   * \sa \ref mcMd_analyzer_MdEnergyAnalyzer_page "param file format"
   *
   * \ingroup McMd_MdAnalyzer_Module
   */
   class MdEnergyAnalyzer : public AverageListAnalyzer<MdSystem>
   {

   public:
   
      /**
      * Constructor.
      *
      * \param system  parent MdSystem object. 
      */
      MdEnergyAnalyzer(MdSystem& system);
   
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

   protected:

      /**
      * Compute values of energy components, store in values_ array.
      */  
      void compute();
  
 
   };

}
#endif 
