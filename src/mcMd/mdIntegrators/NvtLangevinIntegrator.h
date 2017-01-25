#ifndef MCMD_NVT_LANGEVIN_INTEGRATOR_H
#define MCMD_NVT_LANGEVIN_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * A NVT molecular dynamics integrator with a Langevin thermostat.
   *
   * \sa \ref mcMd_integrator_NvtLangevinIntegrator_page "parameter file format"
   * \sa \ref algorithms_Langevin_page "algorithm"
   *
   * This class approximately integrates the Langevin equation
   * \f[
   *    m\frac{d{\bf v}}{dt}  
   *    = -\frac{\partial U}{\partial {\bf r}}
   *   - m\gamma {\bf v} + {\bf f}^{\rm (r)} ,
   * \f]
   * in which \f$\gamma\f$ is a velocity relaxation rate (inverse 
   * time) parameter and \f${\bf f}^{\rm (r)}\f$ is a random force.
   * 
   * \ingroup McMd_MdIntegrator_Module
   */
   class NvtLangevinIntegrator : public MdIntegrator
   {
   
   public:

      /**
      * Constructor. 
      *
      * \param system parent MdSystem
      */
      NvtLangevinIntegrator(MdSystem& system);

      /**
      * Destructor.   
      */
      virtual ~NvtLangevinIntegrator();

      /**
      * Read parameters from file and initialize this MdSystem.
      *
      * \param in input file stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load the internal state to an archive.
      *
      * \param ar archive object.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save the internal state to an archive.
      *
      * \param ar archive object.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Setup private variables before main loop.
      */
      virtual void setup();

      /**
      * Take a complete NVE MD integration step.
      */
      virtual void step();

   private:

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> prefactors_;

      /// Constant for friction force.
      DArray<double> cv_;

      /// Constant for random force.
      DArray<double> cr_;

      /// Velocity autocorrelation decay rate.
      double gamma_;

   };

} 
#endif
