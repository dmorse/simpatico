#ifndef ANGLE_FACTORY_H
#define ANGLE_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class System;
   class AnglePotential;

   /**
   * Factory for subclasses MdAnglePotential or McAnglePotential.
   * 
   * \ingroup Angle_Module
   */
   class AngleFactory : public Factory<AnglePotential>
   {
   
   public:
   
      /**
      * Default constructor.
      */
      AngleFactory(System& system);

      /**
      * Return a pointer to a new AnglePotential, if possible.
      */
      AnglePotential* factory(const std::string& subclass) const;

   private:

      // Pointer to parent system.
      System* systemPtr_;

   };
  
}

#endif // ifndef ANGLE_FACTORY_H
