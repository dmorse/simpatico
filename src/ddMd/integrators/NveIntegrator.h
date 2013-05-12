#ifndef DDMD_NVE_INTEGRATOR_H
#define DDMD_NVE_INTEGRATOR_H

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
   class NveIntegrator : public TwoStepIntegrator
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

      // Time step
      double  dt_;
  
      // Factors of 0.5*dt_/mass, calculated in setup() 
      DArray<double> prefactors_;      

   };

}
#endif
