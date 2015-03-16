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

   class System;
   using namespace Util;

   class LinearGenerator : public Generator
   {

   public:

      LinearGenerator(System& system);

      virtual void generate(int speciesId, int nMolecule,
                            DArray<double> exclusionRadius);

   private:

      bool placeAtom(Molecule& molecule, int atomId,
                             DArray<double> exclusionRadius);
      
      bool placeMolecule(Molecule& molecule, DArray<double> exclusionRadius);

   };

}
#endif
