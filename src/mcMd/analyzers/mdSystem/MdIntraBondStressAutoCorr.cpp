/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/mdSystem/MdIntraBondStressAutoCorr.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/analyzers/system/IntraBondStressAutoCorr.tpp>

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
