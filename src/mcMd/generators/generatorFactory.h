#ifndef MCMD_GENERATOR_FACTORY_H
#define MCMD_GENERATOR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Generator;
   class Species;
   class System;

   Generator* generatorFactory(Species& species, System& system);

}
#endif
