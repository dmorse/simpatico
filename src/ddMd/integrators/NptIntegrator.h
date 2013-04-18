#ifndef DDMD_NPT_INTEGRATOR_H
#define DDMD_NPT_INTEGRATOR_H

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
   * A reversible symplectic NPT integrator.
   *
   * Algorithm is based on Martyna,Tobias and Klein (J. Chem. Phys. 101, 4177-4189, 1994).
   *
   * \ingroup DdMd_Integrator_Module
   */
   class NptIntegrator : public TwoStepIntegrator
   {

   public:

      /**
      * Constructor.
      */
      NptIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~NptIntegrator();

      /**
      * Read required parameters.
      *
      * For velocity-verlet algorithm, reads the time step dt.
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
      * Initialize internal state variables xi, eta, and nu to zero.
      */
      virtual void initDynamicalState();

      /**
      * Setup state just before main loop.
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

   private:

      double xi_;
      double eta_;
      Vector nu_;
      double  dt_;
      double tauT_;
      double tauP_;
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
