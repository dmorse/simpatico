#ifndef DDMD_NVE_INTEGRATOR_H
#define DDMD_NVE_INTEGRATOR_H

#include "Integrator.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;
   using namespace Util;

   /**
   * A velocity-Verlet constant energy integrator.
   *
   * All operations of this class are local (no MPI). 
   */
   class NveIntegrator : public Integrator
   {

   public:

      /**
      * Constructor.
      */
      NveIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~NveIntegrator();

      /**
      * Read required parameters.
      *
      * For velocity-verlet algorithm, reads the time step dt.
      */
      void readParam(std::istream& in);

      /**
      * Setup state just before integration.
      */
      void setup();

      /**
      * Implement one MD step.
      */
      void step();

   private:

      double  dt_;
   
      DArray<double> prefactors_;      

   };

}
#endif
