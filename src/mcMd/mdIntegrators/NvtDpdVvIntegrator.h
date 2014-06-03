#ifndef MCMD_NVT_DPD_VV_INTEGRATOR_H
#define MCMD_NVT_DPD_VV_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <util/containers/DArray.h>
#include <util/space/Vector.h>

#include <iostream>

namespace Util
{
   class Random;
   class EnergyEnsemble;
}

namespace McMd
{

   using namespace Util;

   class PairList;

   /**
   * A velocity-Verlet dissipative particle dynamics (DPD) integrator.
   *
   * This class implements a simple velocity-Verlet (VV) algorithm for
   * the dissipitative particle dynamics (DPD) equations of motion.
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class NvtDpdVvIntegrator : public MdIntegrator
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent MdSystem object
      */
      NvtDpdVvIntegrator(MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~NvtDpdVvIntegrator();

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
      * Take a complete integration step.
      */
      virtual void step();

   private:

      static const int MaxAtomType;

      /// Dissipative forces for atoms.
      DArray<Vector> dissipativeForces_;

      /// Dissipative forces for atoms.
      DArray<Vector> randomForces_;

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> dtMinvFactors_;

      /// Cutoff for DPD thermostat.
      double cutoff_;

      /// DPD pair drag coefficient.
      double gamma_;

      /// Prefactor for random forces.
      double sigma_;

      /// Target absolute temperature (units of energy).
      double temperature_;

      /// Square of cutoff_.
      double cutoffSq_;

      /// Pointer to pair list
      const PairList* pairListPtr_;

      /// Pointer to boundary object.
      const Boundary* boundaryPtr_;

      /// Pointer to random object.
      Random* randomPtr_;

      /// Pointer to random object.
      const EnergyEnsemble* energyEnsemblePtr_;

      /// Atom capacity of parent simulation.
      int atomCapacity_;

      /// Has this object been initialized?
      bool isInitialized_;

      /*
      * Calculate random and dissipative DPD forces.
      *
      * This function computes and stores dissipative forces on all
      * atoms, and computes and random forces iff computeRandom is
      * true. Dissipative and random forces are stored in separate
      * arrays from the conservative atomic forces.
      *
      * \param computeRandom If true compute new random forces.
      */
      void computeDpdForces(bool computeRandom);

   };

}
#endif
