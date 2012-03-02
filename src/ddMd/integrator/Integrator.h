#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <util/param/ParamComposite.h>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class System;
   using namespace Util;

   /**
   * An Integrator numerically integrates the equations of motion. 
   * 
   * This class implements the velocity-Verlet algorithm. The 
   * integrator uses the public interface of a parent System to
   * evaluated forces and access atomic forces and velocities.
   *
   * All operations of this class are local (no MPI). 
   */
   class Integrator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Integrator(System& system);

      /**
      * Destructor.
      */
      ~Integrator();

      /**
      * Read required parameters.
      *
      * For velocity-verlet algorithm, reads the time step dt.
      */
      void readParam(std::istream& in);

      /**
      * Initialize just before integration.
      */
      void initialize();

      /**
      * Implement one step.
      */
      void step();

   private:

      System* systemPtr_;

      double  dt_;

   };

}
#endif
