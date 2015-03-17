#ifndef MCMD_POINT_GENERATOR_H
#define MCMD_POINT_GENERATOR_H

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

   class PointGenerator : public Generator
   {

   public:

      PointGenerator(Species& species, System& system);

   protected:

      bool attemptPlaceMolecule(Molecule& molecule, 
                                const DArray<double>& diameters,
                                CellList& cellList);

   };

}
#endif
