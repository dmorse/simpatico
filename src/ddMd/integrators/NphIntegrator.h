#ifndef DDMD_NPH_INTEGRATOR_H
#define DDMD_NPH_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TwoStepIntegrator.h"      // base class

namespace DdMd
{

   class Simulation;
   using namespace Util;

   /**
   * A velocity-Verlet constant energy integrator.
   *
   * \ingroup DdMd_Integrator_Module
   */
   class NphIntegrator : public TwoStepIntegrator
   {

   public:

      /**
      * Constructor.
      */
      NphIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~NphIntegrator();

      /**
      * Read required parameters.
      *
      * For velocity-verlet algorithm, reads the time step dt.
      */
      void readParam(std::istream& in);

      /**
      * Setup state just before integration.
      *
      * Must be called before run() method.
      */
      void setup();

   protected:

      /**
      * Execute first step of two-step integrator.
      */
      virtual void integrateStep1();

      /**
      * Execute secodn step of two-step integrator.
      */
      virtual void integrateStep2();

   private:
      Vector nu_;
      double  dt_;
      double W_;
      double T_target_;
      double P_target_;
      double T_kinetic_;
      Vector P_curr_diag_;
      LatticeSystem mode_;
      unsigned int ndof_;

      Vector exp_v_fac_;
      Vector sinhx_fac_v_;
      double mtk_term_2_;
      DArray<double> prefactors_;

   };

}
#endif
