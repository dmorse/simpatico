#ifndef DDMD_INTEGRATOR_H
#define DDMD_INTEGRATOR_H

#include <util/param/ParamComposite.h>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * An Integrator numerically integrates the equations of motion. 
   *
   * \ingroup DdMd_Integrator_Module
   */
   class Integrator : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Integrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~Integrator();

      #if 0
      /**
      * Read required parameters.
      *
      * For velocity-verlet algorithm, reads the time step dt.
      */
      virtual void readParam(std::istream& in);
      #endif

      /**
      * Initialize just before integration.
      */
      virtual void setup() = 0;

      /**
      * Implement one step.
      */
      virtual void step() = 0;

   protected:

      /**
      * Get reference to parent Simulation.
      */ 
      Simulation& simulation();

   private:

      Simulation* simulationPtr_;

   };

   /// Get reference to parent Simulation.
   inline Simulation& Integrator::simulation() 
   {  return *simulationPtr_; }

}
#endif
