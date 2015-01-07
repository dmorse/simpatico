#ifndef DDMD_NVT_LANGEVIN_INTEGRATOR_H
#define DDMD_NVT_LANGEVIN_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
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
   * The time-stepping is a conventional two-step velocity-Verlet 
   * integrator similar to that of NveVvIntegrator, with a total
   * force of the form:
   * \f{eqnarray*}
   *    {\bf f} = -\frac{\partial U}{\partial {\bf r}}
   *            + c_{v}c_{v}{\bf v} + c_{r}{\bf u} .
   * \f}
   * in \f$c_{v} < 0\f$ and \f$c_{r}\f$ are numerical coefficients that
   * must be chosen to produce target temperature, as discussed below,
   * and \f${\bf u}\f$ is a dimensioness random number whose Cartesian 
   * components are each uniformly distributed random numbers in the range 
   * \f$[-1/2, 1/2]\f$.  The drag and random forces are both evaluated and 
   * added to the total force at the beginning of the second step of the 
   * two-step velocity-Verlet integration, immediately after evaluation of 
   * all conservative forces. The velocity used to calculate this force is 
   * thus the "half-updated" velocity produced by the first half of the 
   * two-step velocity-Verlet integrator.
   * 
   * Values of coefficients:
   * In the absence of any deterministic force, the velocity change between 
   * force evaluations is given by
   * \f[
   *    {\bf v} \rightarrow ( 1 + c_{v} \Delta t/m ){\bf v} +
   *             (c_{r}\Delta t/m) {\bf u} 
   * \f]
   * By convention, we choose the constant \f$c_{v}\f$ such that 
   * \f$1 + c_{v}\Delta t/m = e^{-\gamma\Delta t}\f$, where \f$\gamma\f$ is 
   * the decay rate parameter input by the user. This guarantees that the
   * decay of velocity correlations between subsequent steps is that given 
   * by an exact solution of the Langevin equation. This gives:
   * \f[
   *     c_{v} = -m (1 - e^{-\gamma\Delta t})/\Delta t .
   * \f]
   * The coefficient \f$c_{r}\f$ for the random force is then chosen so 
   * that, in the absence of conservative forces, the mapping would exactly 
   * preserve a variance 
   * \f[ 
   *     \langle v_{i} v_{j} \rangle = \delta_{ij} d k_{B}T/m . 
   * \f] 
   * for Cartesian components of the midstep velocities at temperature 
   * \f$T\f$.  The coefficient \f$d\f$ is a dimensionless constant that 
   * approaches unity in the limit \f$\gamma\Delta t \rightarrow 0\f$, 
   * but that is chosen here such that the velocity at the beginning 
   * and end of each step (rather than at the midstep) satisfies the 
   * usual equipartition theorem, with d=1. Applying this criteria, and
   * using the identity \f$\langle u_{i} u_{j} \rangle = \delta_{ij}/12\f$ 
   * for the Cartesian components of \f${\bf u}\f$, we obtain:
   * \f[
   *     c_{r} = \sqrt{12 m k_{B}T d ( 1 - e^{-2\gamma\Delta t})}/\Delta t
   * \f]
   * Requiring that the end-of-step velocity satisfy equipartitions in 
   * the absence of external forces then yields a coefficient
   * \f[
   *    d = 2/(1 + e^{-\gamma\Delta t}) .
   * \f]
   * 
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
