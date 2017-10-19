/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ColVarPotentialTmpl.h"
#include <mcMd/simulation/System.h> 

namespace McMd
{

   using namespace Util;

   /**
   * Constructor.
   */
   template <class ColVarType, class BiasType>
   ColVarPotentialTmpl<ColVarType, BiasType>::ColVarPotentialTmpl(System& system)
    : SpecialPotential(false),
      systemPtr_(&system),
      colVar_(),
      bias_(),
   {  setClassName("ColVarPotential"); }

   /*
   * Destructor.
   */
   template <class ColVarType, class BiasType>
   ColVarPotentialTmpl<ColVarType, BiasType>::~ColVarPotentialTmpl()
   {}

   /*
   * Read parameters from file.
   */
   template <class ColVarType, class BiasType>
   void 
   ColVarPotentialTmpl<ColVarType, BiasType>::readParameters(std::istream& in)
   {
      colVar_.readParam(in);
      bias_.readParam(in);
   }

   /*
   * Compute total energy.
   */
   template <class ColVarType, class BiasType>
   void ColVarPotentialTmpl<ColVarType, BiasType>::computeEnergy()
   {
      double cv = colVar_.value();
      return bias_.value(cv);
   }

   /*
   * Unset energy and value of collective variable.
   */
   template <class ColVarType, class BiasType> void 
   ColVarPotentialTmpl<ColVarType, BiasType>::unsetEnergy()
   {
      SpecialPotential::unsetEnergy();
      colVar_.unset();
   }

   /*
   * Add forces from this potential to all atomic forces.
   */
   template <class ColVarType, class BiasType>
   void ColVarPotentialTmpl<ColVarType, BiasType>::addForces()
   {  
      double cv = colVar_.value();
      double dwdc = bias_.derivative(cv);
      colVar.addForces(dwdc);
   }

} 
