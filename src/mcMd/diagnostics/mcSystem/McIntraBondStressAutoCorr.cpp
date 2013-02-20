#ifndef MCMD_MC_INTRA_BOND_STRESS_AUTO_CORR_CPP
#define MCMD_MC_INTRA_BOND_STRESS_AUTO_CORR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/mcSystem/McIntraBondStressAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/diagnostics/system/IntraBondStressAutoCorr.tpp>

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
#endif

