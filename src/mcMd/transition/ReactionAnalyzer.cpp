#ifndef MCMD_REACTION_ANALYZER_H
#define MCMD_REACTION_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   ReactionAnalyzer::ReactionAnalyzer(MdSystem& system);

   void ReactionAnalyzer::readParameter(std::istream& in)
   {
      bool checkBase = false;
      readInterval(in, checkBase);
      read<double>(in, "reactantCoordinate", reactantCoordinate_);
      read<double>(in, "productCoordinate", productCoordinate_);
   }

   void ReactionAnalyzer::sample(int iStep)
   {
      currentCoordinate_ = computeReactionCoordinate();
      iStep_ = iStep;
   }

}
#endif
