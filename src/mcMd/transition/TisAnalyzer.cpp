#ifndef MCMD_TIS_ANALYZER_H
#define MCMD_TIS_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   TisReactionCoordinate(MdSystem& system);

   void readParameter(std::istream& in)
   {
   }

   void setLower(double lowerRC);
   {  lowerRC_ = lowerRC_; }

   void setUpper(double upperRC);
   {  upperRC_ = upperRC_; }

}
#endif
