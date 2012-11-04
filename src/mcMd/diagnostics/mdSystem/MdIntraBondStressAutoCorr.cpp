#ifndef MCMD_MD_INTRA_BOND_STRESS_AUTO_CORR_CPP
#define MCMD_MD_INTRA_BOND_STRESS_AUTO_CORR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/mdSystem/MdIntraBondStressAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/diagnostics/system/IntraBondStressAutoCorr.tpp>

namespace McMd
{

   using namespace Util;
  
   /*
   * Constructor.
   */
   MdIntraBondStressAutoCorr::MdIntraBondStressAutoCorr(MdSystem &system)
    : IntraBondStressAutoCorr<MdSystem>(system)
   {  setClassName("MdIntraBondStressAutoCorr"); }

   /*
   * Destructor.
   */
   MdIntraBondStressAutoCorr::~MdIntraBondStressAutoCorr()
   {}

}
#endif

