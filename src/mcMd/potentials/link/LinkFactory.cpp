#ifdef  MCMD_LINK
#ifndef LINK_FACTORY_CPP
#define LINK_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/link/LinkFactory.h>
#include <mcMd/simulation/System.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/potentials/link/LinkPotentialImpl.h>

// Bond Potential evaluator classes
#include <mcMd/potentials/bond/HarmonicBond.h>
#include <mcMd/potentials/bond/HarmonicL0Bond.h>
#include <mcMd/potentials/bond/FeneBond.h>

namespace McMd
{

   /*
   * Default constructor.
   */
   LinkFactory::LinkFactory(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a new BondPotential, if possible.
   */
   BondPotential* LinkFactory::factory(const std::string& name) const
   {
      BondPotential* ptr = 0;
      if (name == "HarmonicBond") {
         ptr = new LinkPotentialImpl<HarmonicBond>(*systemPtr_);
      } else
      if (name == "HarmonicL0Bond") {
         ptr = new LinkPotentialImpl<HarmonicL0Bond>(*systemPtr_);
      } else
      if (name == "FeneBond") {
         ptr = new LinkPotentialImpl<FeneBond>(*systemPtr_);
      }

      return ptr;
   }

}
#endif
#endif
