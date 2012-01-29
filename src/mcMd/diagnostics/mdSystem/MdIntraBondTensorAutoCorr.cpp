#ifndef MD_INTRA_BOND_TENSOR_AUTO_CORR_CPP
#define MD_INTRA_BOND_TENSOR_AUTO_CORR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/mdSystem/MdIntraBondTensorAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/diagnostics/system/IntraBondTensorAutoCorr_inc.h>

namespace McMd
{

   using namespace Util;
  
   /*
   * Constructor.
   */
   MdIntraBondTensorAutoCorr::MdIntraBondTensorAutoCorr(MdSystem &system)
    : IntraBondTensorAutoCorr<MdSystem>(system)
   {}

   /*
   * Destructor.
   */
   MdIntraBondTensorAutoCorr::~MdIntraBondTensorAutoCorr()
   {}

}
#endif

