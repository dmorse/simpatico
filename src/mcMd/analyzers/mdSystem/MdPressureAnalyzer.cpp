/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPressureAnalyzer.h"
#include <mcMd/analyzers/base/AverageListAnalyzer.tpp>
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>

#include <util/accumulators/Average.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdPressureAnalyzer::MdPressureAnalyzer(MdSystem& system) 
    : AverageListAnalyzer<MdSystem>(system),
      virialId_(-1),
      kineticId_(-1),
      totalId_(-1)
   {  setClassName("MdPressureAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdPressureAnalyzer::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<MdSystem>::readParameters(in);

      // Set id values.
      virialId_ = 0;
      kineticId_ = 1;
      totalId_ = 2;
      initializeAccumulators(3);

      // Set names
      setName(virialId_, "virial");
      setName(kineticId_, "kinetic");
      setName(totalId_, "total");
   }

   /*
   * Load parameters from archive when restarting. 
   */
   void MdPressureAnalyzer::loadParameters(Serializable::IArchive& ar) 
   {
      AverageListAnalyzer<MdSystem>::loadParameters(ar);
      virialId_ = 0;
      kineticId_ = 1;
      totalId_ = 2;
   }

   /*
   * Output energy to file
   */
   void MdPressureAnalyzer::compute() 
   {
      MdSystem& sys = system();

      double virial;
      sys.computeVirialStress<double>(virial);
      setValue(virialId_, virial);

      double kinetic; 
      sys.computeKineticStress<double>(kinetic);
      setValue(kineticId_, kinetic);

      double total = kinetic + virial;
      setValue(totalId_, total);
   }

}
