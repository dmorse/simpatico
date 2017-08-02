/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdExternalFacade.h"

namespace McMd
{

   using namespace Util;

   /**
   * Constructor.
   */
   MdExternalFacade::MdExternalFacade(System& system)
    : MdPotentialFacade<ExternalPotential, ExternalFactory>(system)
   {  setClassName("MdExternalFacade"); }

} 
