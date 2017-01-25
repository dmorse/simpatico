/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/mdSystem/MdIntraBondTensorAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/analyzers/system/IntraBondTensorAutoCorr.tpp>

namespace McMd
{

   using namespace Util;
  
   /*
   * Constructor.
   */
   MdIntraBondTensorAutoCorr::MdIntraBondTensorAutoCorr(MdSystem &system)
    : IntraBondTensorAutoCorr<MdSystem>(system)
   { setClassName("MdIntraBondTensorAutoCorr"); }

   /*
   * Destructor.
   */
   MdIntraBondTensorAutoCorr::~MdIntraBondTensorAutoCorr()
   {}

}
