/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairEnergyAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/potentials/pair/PairPotential.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   PairEnergyAnalyzer::PairEnergyAnalyzer(Simulation& simulation) 
    : AverageAnalyzer(simulation)
   {  setClassName("PairEnergyAnalyzer"); }

   /*
   * Destructor.
   */
   PairEnergyAnalyzer::~PairEnergyAnalyzer() 
   {}  

   /*
   * Read interval and outputFileName. 
   */
   void PairEnergyAnalyzer::readParameters(std::istream& in) 
   {
      AverageAnalyzer::readParameters(in);
      readFArray<int, 2>(in, "typeIdPair", typeIdPair_);
   }

   /*
   * Load internal state from an archive.
   */
   void PairEnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      AverageAnalyzer::loadParameters(ar);
      loadFArray<int, 2>(ar, "typeIdPair", typeIdPair_);
   }

   /*
   * Load internal state from an archive.
   */
   void PairEnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      AverageAnalyzer::save(ar);
      ar << typeIdPair_;
   }

   /*
   * Compute current value.
   */
   void PairEnergyAnalyzer::compute() 
   {  
      //simulation().computePairEnergies(); 
      MPI::Intracomm& communicator = simulation().domain().communicator();  
      simulation().pairPotential().computePairEnergies(communicator); 
   }

   double PairEnergyAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      int i = typeIdPair_[0];
      int j = typeIdPair_[1];
      DMatrix<double> pair = simulation().pairPotential().pairEnergies();
      return 0.5*(pair(i,j) + pair(j,i));
   }

}
