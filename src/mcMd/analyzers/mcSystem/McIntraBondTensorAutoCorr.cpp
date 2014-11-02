/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/mcSystem/McIntraBondTensorAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/analyzers/system/IntraBondTensorAutoCorr.tpp>

namespace McMd
{

   using namespace Util;
  
   /*
   * Constructor.
   */
   McIntraBondTensorAutoCorr::McIntraBondTensorAutoCorr(McSystem &system)
    : IntraBondTensorAutoCorr<McSystem>(system)
   {  setClassName("McIntraBondTensorAutoCorr");  }

   /*
   * Destructor.
   */
   McIntraBondTensorAutoCorr::~McIntraBondTensorAutoCorr()
   {}

}
