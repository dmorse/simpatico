#ifdef  MCMD_ANGLE
#ifndef ANGLE_FACTORY_CPP
#define ANGLE_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/angle/AngleFactory.h>

// AnglePotential base and implementation classes
#include <mcMd/potentials/angle/AnglePotential.h>
#include <mcMd/potentials/angle/AnglePotentialImpl.h>

// Angle Potential evaluator classes
#include <mcMd/potentials/angle/CosineAngle.h>
#include <mcMd/potentials/angle/CosineSqAngle.h>

namespace McMd
{

   /**
   * Default constructor.
   */
   AngleFactory::AngleFactory(System& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a new AnglePotential, if possible.
   */
   AnglePotential* 
   AngleFactory::factory(const std::string& name) const
   {
      AnglePotential* ptr = 0;

      ptr = trySubfactories(name);
      if (ptr) return ptr;

      if (name == "CosineAngle") {
         ptr = new AnglePotentialImpl<CosineAngle>(*systemPtr_);
      } else
      if (name == "CosineSqAngle") {
         ptr = new AnglePotentialImpl<CosineSqAngle>(*systemPtr_);
      }
      return ptr;
   }

}
#endif // ifndef ANGLE_FACTORY_CPP
#endif // ifdef MCMD_ANGLE
