#ifndef DDMD_NVT_INTEGRATOR_H
#define DDMD_NVT_INTEGRATOR_H

#include "TwoStepIntegrator.h"

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   class Simulation;
   using namespace Util;

   /**
   * A Nose-Hoover constant temperature, constant volume integrator.
   *
   * \ingroup DdMd_Integrator_Module
   */
   class NvtIntegrator : public TwoStepIntegrator
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
      */
      void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

   protected:

      /**
      * Setup state just before integration.
      */
      void setup();

      /**
      * Execute first step of two-step integrator.
      */
      virtual void integrateStep1();

      /**
      * Execute secodn step of two-step integrator.
      */
      virtual void integrateStep2();

      /**
      * Initialize internal dynamical state variables to default value.
      */
      virtual void initDynamicalState();

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
