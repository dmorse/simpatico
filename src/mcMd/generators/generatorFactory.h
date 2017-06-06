#ifndef MCMD_GENERATOR_FACTORY_H
#define MCMD_GENERATOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Simp
{
   class Species;
}

namespace McMd
{

   class Generator;
   class McSystem;
   class MdSystem;

   using namespace Simp;

   /**
   * Instantiates generator for on species in an McSystem.
   *
   * Returns a pointer to the new generator, or null if no appropriate
   * generator is found.
   *
   * \ingroup McMd_Generator_Module
   * 
   * \param species Species object
   * \param system parent McSystem
   */
   Generator* generatorFactory(Species& species, McSystem& system);

   /**
   * Instantiates generator for on species in an MdSystem.
   *
   * Returns a pointer to the new generator, or null if no appropriate
   * generator is found.
   * 
   * \ingroup McMd_Generator_Module
   * 
   * \param species Species object
   * \param system parent McSystem
   */
   Generator* generatorFactory(Species& species, MdSystem& system);

}
#endif
