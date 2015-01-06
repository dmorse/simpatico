#ifndef MCMD_NVT_LANGEVIN_INTEGRATOR_H
#define MCMD_NVT_LANGEVIN_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <iostream>

namespace Util{ class EnergyEnsemble; } 

namespace McMd
{

   using namespace Util;

   /**
   * A Langevin NVT molecular dynamics integrator.
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class NvtLangevinIntegrator : public MdIntegrator
   {
   
   public:

      /// Constructor. 
      NvtLangevinIntegrator(MdSystem& system);
 
      /// Destructor.   
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
      DArray<double> vcoeff_;

      /// Constant for random force.
      DArray<double> fcoeff_;

      /// Velocity autocorrelation decay rate.
      double gamma_;

   };

} 
#endif
