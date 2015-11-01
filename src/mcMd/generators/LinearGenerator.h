#ifndef MCMD_LINEAR_GENERATOR_H
#define MCMD_LINEAR_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Generator.h"

namespace McMd
{

   class Species;
   class System;
   class Molecule;
   class CellList;
   using namespace Util;

   /**
   * Generates random configurations for linear molecules.
   *
   * \ingroup McMd_Generators_Module
   *
   */
   class LinearGenerator : public Generator
   {

   public:

      LinearGenerator(Species& species, System& system);

   protected:

      bool attemptPlaceMolecule(Molecule& molecule, 
                                const DArray<double>& diameters,
                                CellList& cellList);

   };

}
#endif
