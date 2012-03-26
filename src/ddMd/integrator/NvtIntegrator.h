#ifndef DDMD_NVT_INTEGRATOR_H
#define DDMD_NVT_INTEGRATOR_H

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
   * \ingroup DdMd_Integrator_Module
   */
   class NvtIntegrator : public Integrator
   {

   public:

      /**
      * Constructor.
      */
      NvtIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~NvtIntegrator();

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

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> prefactors_;

      /// Time step.
      double  dt_;
   
      /// Target temperature
      double T_target_;

      /// Current temperature from kinetic energy
      double T_kinetic_;

      /// Nose-Hover thermostat scaling variable.
      double xi_;

      /// Time derivative of xi
      double xiDot_;

      /// Relaxation time for energy fluctuations.
      double tauT_;

      /// Relaxation rate for energy fluctuations.
      double nuT_;

      /// Total number of atoms in simulation.
      int nAtom_;

   };

}
#endif
