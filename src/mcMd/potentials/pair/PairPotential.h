#ifndef MCMD_PAIR_POTENTIAL_H
#define MCMD_PAIR_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/misc/EnergyCalculator.h>
#include <mcMd/potentials/misc/StressCalculator.h>

// #include <util/space/Tensor.h>
// #include <util/misc/Setable.h>

#include <util/global.h>
#include <string>

namespace Util
{
   class Vector;
}

namespace McMd
{

   using namespace Util;

   /**
   * Interface for a Pair Potential.
   *
   * \ingroup McMd_Pair_Module
   */
   class PairPotential : public EnergyCalculator, public StressCalculator
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~PairPotential()
      {}

      using EnergyCalculator::energy;

      /// \name Pair Interaction Interface
      //@{ 

      /**
      * Return pair energy for a single pair.
      */
      virtual 
      double energy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double forceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  new value of parameter
      */
      virtual void set(std::string name, int i, int j, double value) = 0;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      virtual double get(std::string name, int i, int j) const = 0;

      /**
      * Return maximum cutoff distance.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Return name of pair interaction class (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const = 0;

      //@}

   };

}
#endif
