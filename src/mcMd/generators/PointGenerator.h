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

   /**
   * Generator for monoatomic molecules (atoms).
   *
   * \ingroup McMd_Generator_Module
   */
   class PointGenerator : public Generator
   {

   public:

      /**
      * Constructor.
      *
      * \param species  molecular Species to be generated
      * \param system  parent System object
      */
      PointGenerator(Species& species, System& system);

   protected:

      /**
      * Attempt to place an "molecule" (i.e., an atom).
      *
      * If successful, the atom is added to the CellList.
      *
      * \param molecule reference to Molecule object
      * \param diameters  array of hard-core exclusion diameters
      * \param cellList  CellList object, modified if successful.
      * \return true for success, false for failure
      */
      bool attemptPlaceMolecule(Molecule& molecule, 
                                const DArray<double>& diameters,
                                CellList& cellList);

   };

}
#endif
