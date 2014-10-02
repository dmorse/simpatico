#ifndef DDMD_ANGLE_FACTORY_H
#define DDMD_ANGLE_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>                  // base class template
#include <ddMd/potentials/angle/AnglePotential.h>  // template argument

#include <string>

namespace DdMd
{

   class Simulation;

   /**
   * Factory for AnglePotential objects.
   *
   * \ingroup DdMd_Angle_Module
   */
   class AngleFactory : public Factory<AnglePotential>
   {

   public:
   
      /**
      * Default constructor.
      */
      AngleFactory(Simulation& simulation);

      /**
      * Return a pointer to a new McAngleInteration, if possible.
      */
      AnglePotential* factory(const std::string& subclass) const;

   private:

      // Pointer to the parent Simulation.
      Simulation* simulationPtr_;

   };
  
}
#endif
