#ifndef MCMD_NPH_INTEGRATOR_H
#define MCMD_NPH_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * An explictly reversible/measure-preserving Parrinello-Rahman type NPH integrator.
   *
   * TwoStepNPH can be used in three modes:
   *  - isotropic volume flucutations, cubic box
   *  - anisotropic volume fluctuations, orthorhombic box
   *    (three independent lengths)
   *  - anisotropic volume fluctuations, tetragonal box
   *    (two independent lengths, Ly = Lz)
   *
   *  The first mode is an implementation of the Andersen barostat,
   *  the second and the third mode are based on the adapted Parrinello-Rahman
   *  equations of motion.
   *
   *  The integration scheme used to implement the equations of motion is
   *  explicitly reversible and measure-preserving. It is based
   *  on the Trotter expansion technique introduced by Tuckerman et al. J. Chem.
   *  Phys.  97, pp. 1990 (1992).
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class NphIntegrator : public MdIntegrator
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent MdSystem object
      */
      NphIntegrator(MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~NphIntegrator();

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

      /**
      * Get the barostat energy.
      */
      virtual double barostatEnergy();

      /**
      * Get the barostat mass.
      */
      virtual double barostatMass() const;

      /**
      * Get the integrator mode.
      */
      virtual LatticeSystem mode() const;

      /**
      * Get the barostat momentum.
      */
      //virtual Vector& eta();

      /**
      * Set the barostat momentum component of an index to a value eta
      *
      * \param index index of momentum component.
      * \param eta value to which momentum component must be set.
      */
      virtual void setEta(unsigned int index, double eta);

   private:

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> prefactors_;

      /// internal state variables (momentum conjugate to the box lengths/volume)
      Vector eta_;

      /// barostat mass
      double W_;

      /// current stress tensor (diagonal elements)
      Vector currP_;

      /// integration mode
      LatticeSystem mode_;

   };

}
#endif
