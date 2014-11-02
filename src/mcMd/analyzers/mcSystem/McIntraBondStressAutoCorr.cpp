/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/mcSystem/McIntraBondStressAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/analyzers/system/IntraBondStressAutoCorr.tpp>

namespace McMd
{

   using namespace Util;
  
   /*
   * Constructor.
   */
   McIntraBondStressAutoCorr::McIntraBondStressAutoCorr(McSystem &system)
    : IntraBondStressAutoCorr<McSystem>(system)
   {  setClassName("McIntraBondStressAutoCorr"); }

   /*
   * Destructor.
   */
   McIntraBondStressAutoCorr::~McIntraBondStressAutoCorr()
   {}

}
