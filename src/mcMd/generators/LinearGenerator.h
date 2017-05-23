#ifndef MCMD_LINEAR_GENERATOR_H
#define MCMD_LINEAR_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Generator.h"

namespace Simp
{
   class Species;
}

namespace McMd
{

   class System;
   class Molecule;
   class CellList;

   using namespace Util;
   using namespace Simp;

   /**
   * Generates random configurations for linear molecules.
   *
   * \ingroup McMd_Generator_Module
   *
   */
   class LinearGenerator : public Generator
   {

   public:

      /**
      * Constructor.
      *
      * \param species molecular Species to be generated
      * \param system  parent System object
      */
      LinearGenerator(Species& species, System& system);

   protected:

      /**
      * Attempt to place an entire linear chain.
      *
      * If successful, all atoms are added to the CellList.
      * If unsuccesful, CellList is unchanged.
      *
      * \param molecule reference to Molecule object
      * \param diameters  array of hard-core exclusion diameters
      * \param cellList  CellList object, modified if successful.
      * \return true for success, false for failure
      */
      bool attemptPlaceMolecule(Molecule& molecule, 
                                Array<double> const & diameters,
                                CellList& cellList);

   };

}
#endif
