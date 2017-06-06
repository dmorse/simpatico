/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/transition/ReactionAnalyzer.h>

namespace McMd 
{

   // using namespace Util;

   ReactionAnalyzer::ReactionAnalyzer(MdSystem& system)
    : SystemAnalyzer<MdSystem>(system)
   {}

   void ReactionAnalyzer::readParameter(std::istream& in)
   {
      bool checkBase = false;
      readInterval(in, checkBase);
      read<double>(in, "reactantCoordinate", reactantCoordinate_);
      read<double>(in, "productCoordinate", productCoordinate_);
   }

   void ReactionAnalyzer::sample(long iStep)
   {
      currentCoordinate_ = computeReactionCoordinate();
      iStep_ = iStep;
   }

}
