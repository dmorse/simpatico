#ifndef MCMD_GENERATOR_H
#define MCMD_GENERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace McMd
{

   using namespace Util;

   class System;
   class McSystem;
   class MdSystem;
   class BondPotential;

   class Generator {

   public:

      Generator(System& system);

      #ifdef INTER_BOND
      void setBondPotential(BondPotential& bondPotential);
      #endif

      void generate(int speciesId, int nMolecule,
                    DArray<double> exclusionRadius);

   protected:

      System& system()
      {  return *systemPtr_; }

      #ifdef INTER_BOND
      BondPotential& bondPotential()
      {  return *bondPotentialPtr_; }
      #endif

   private:

      System* systemPtr_;

      #ifdef INTER_BOND
      BondPotential* bondPotentialPtr_;
      #endif

   };

   class McGenerator : public Generator 
   {
   public:
      McGenerator(McSystem& system);
   };

   class MdGenerator : public Generator 
   {
   public:
      MdGenerator(MdSystem& system);
   };

}
#endif
