#ifndef DDMD_NVT_LANGEVIN_INTEGRATOR_H
#define DDMD_NVT_LANGEVIN_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TwoStepIntegrator.h"      // base class

namespace DdMd
{

   class Simulation;
   using namespace Util;

   /**
   * A NVT molecular dynamics integrator with a Langevin thermostat.
   *
   * This class approximately integrates the Langevin equation
   * \f[
   *    m\frac{d{\bf v}}{dt}  =  
   *    -\frac{\partial U}{\partial {\bf r}}
   *   - m\gamma {\bf v} + {\bf f}^{\rm (r)} ,
   * \f]
   * in which \f$\gamma\f$ is a velocity relaxation rate (inverse 
   * time) parameter, \f${\bf v}\f$ is a particle velocity, and 
   * \f${\bf f}^{\rm (r)}\f$ is a random force.
   * 
   * \sa \ref ddMd_integrator_NvtLangevinIntegrator_page "parameter file format"
   * \sa \ref algorithms_Langevin_page "algorithm"
   * \ingroup DdMd_Integrator_Module
   */
   class NvtLangevinIntegrator : public TwoStepIntegrator
   {

   public:

      /**
      * Constructor.
      */
      NvtLangevinIntegrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~NvtLangevinIntegrator();

      /**
      * Read required parameters.
      *
      * Reads the time step dt.
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
      * Setup state just before main loop.
      *
      * Calls Integrator::setupAtoms(), initializes prefactors_ array.
      */
      void setup();

      /**
      * Execute first step of two-step integrator.
      *
      * Update positions and half-update velocities.
      */
      virtual void integrateStep1();

      /**
      * Execute second step of two-step integrator.
      *
      * Second half-update of velocities.
      */
      virtual void integrateStep2();

   private:

      /// Time step (parameter)
      double  dt_;
  
      /// Velocity autocorrelation decay rate (parameter)
      double gamma_;

      /// Factors of 0.5*dt_/mass, calculated in setup().
      DArray<double> prefactors_;      

      /// Constant for friction force.
      DArray<double> cv_;

      /// Constant for random force.
      DArray<double> cr_;

   };

}
#endif
