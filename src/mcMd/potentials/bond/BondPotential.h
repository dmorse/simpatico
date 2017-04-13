#ifndef MCMD_BOND_POTENTIAL_H
#define MCMD_BOND_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>               // base class
#include <mcMd/potentials/misc/EnergyCalculator.h>   // base class
#include <mcMd/potentials/misc/StressCalculator.h>   // base class

#include <string>

namespace Util
{
   class Vector;
   class Tensor;
   class Random;
}

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * Abstract Bond Potential class.
   *
   * \ingroup McMd_Bond_Module
   */
   class BondPotential : public ParamComposite, 
                         public EnergyCalculator, public StressCalculator
   {

   public:

      /**
      * Constructor .
      */
      BondPotential();

      /**
      * Destructor (does nothing)
      */
      virtual ~BondPotential();

      /// \name Bond Interaction Interface
      //@{

      /**
      * Returns potential energy for one bond.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      */
      virtual double energy(double rSq, int type) const = 0;
   
      // Prevent hiding of inherited function energy();
      using EnergyCalculator::energy;
    
      /**
      * Returns force/distance for one bond, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      * \return     scalar repulsive force divided by distance.
      */
      virtual double forceOverR(double rSq, int type) const = 0;

      /**
      * Return bond length chosen from equilibrium distribution.
      *
      * This function returns a bond length chosen from the Boltzmann
      * distribution of lengths for bonds of random orientation. The
      * distribution P(l) of values of the length l is proportional 
      * to l*l*exp[-beta*phi(l) ], where phi(l) is the bond energy. 
      *
      * \param random pointer to random number generator object.
      * \param beta   inverse temperature
      * \param type   bond type
      * \return bond length chosen from equilibrium distribution.
      */
      virtual 
      double randomBondLength(Util::Random* random, double beta, int type) 
             const = 0;

      /**
      * Modify an interaction parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  bond type index 
      * \param value new value of parameter
      */
      virtual void set(std::string name, int type, double value) = 0;

      /**
      * Get an interaction parameter value, identified by a string.
      *
      * \param name parameter name
      * \param type bond type index 1
      */
      virtual double get(std::string name, int type) const = 0;

      /**
      * Return name of pair interaction class (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name System energy, force, and stress.
      //@{

      /**
      * Calculate the covalent bond energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const = 0;

      /**
      * Add bond forces to all atomic forces.
      */
      virtual void addForces() = 0;

      //@}

   };

} 
#endif
