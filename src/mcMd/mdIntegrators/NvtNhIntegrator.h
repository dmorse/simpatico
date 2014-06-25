#ifndef MCMD_NVT_NH_INTEGRATOR_H
#define MCMD_NVT_NH_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <iostream>

namespace Util{ class EnergyEnsemble; } 

namespace McMd 
{

   using namespace Util;

   /**
   * A Nose-Hoover NVT molecular dynamics integrator.
   *
   * The step() function implements a reversible velocity-verlet MD NVT 
   * integrator step for the Nose-Hoover equations, as described by
   * Winkler, Kraus, and Reineker, J. Chem. Phys. 102, 9018 (1995).
   *
   * Variables names imitate the notation used in the book by
   * D. Frenkel and B. Smit, "Understanding Molecular Simulation,"  
   * Academic Press, 1996. Chapter 6 (Eqs. 6.1.24 - 6.1.27).
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class NvtNhIntegrator : public MdIntegrator
   {
   
   public:

      /**
      * Constructor. 
      *
      * \param system parent MdSystem object
      */
      NvtNhIntegrator(MdSystem& system);

      /** 
      * Destructor.   
      */
      virtual ~NvtNhIntegrator();

      /**
      * Read parameters from file and initialize.
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
      * Setup auxiliary parameters of integrator, just before main loop.
      *
      * This method must be called before the beginning of an MD loop
      * to initialize private variables used by step() method.
      */
      virtual void setup();

      /**
      * Take a complete NVT MD integration step.
      */
      virtual void step();

   private:

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> prefactors_;

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

      /// Pointer to EnergyEnsemble object.
      EnergyEnsemble* energyEnsemblePtr_;

   }; 

} 
#endif
