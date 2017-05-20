#ifdef  SIMP_TETHER
#ifndef MCMD_TETHER_POTENTIAL_H
#define MCMD_TETHER_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
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
   * Abstract Tether Potential class.
   *
   * \ingroup McMd_Tether_Module
   */
   class TetherPotential : public ParamComposite
   {

   public:

      /**
      * Constructor .
      */
      TetherPotential()
      {}

      /**
      * Destructor (does nothing)
      */
      virtual ~TetherPotential()
      {}

      /**
      * Return name of pair evaluator class (e.g., "HarmonicTether").
      */
      virtual std::string evaluatorClassName() const = 0;

      /**
      * Returns potential energy for one tether.
      *
      * \param rSq  square of distance between atom and tether point.
      * \param type type of tether.
      */
      virtual double energy(double rSq, int type) const = 0;
  
      /**
      * Returns force/distance for one tether, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between atom and anchor.
      * \param type type of tether.
      * \return     scalar repulsive force divided by distance.
      */
      virtual double forceOverR(double rSq, int type) const = 0;
  
      /**
      * Calculate the covalent bond energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond potential energy of atom
      */
      virtual double atomEnergy(const Atom& atom) const
      {  UTIL_THROW("Unimplemented method"); }

      /**
      * Add bond forces to all atomic forces.
      */
      virtual void addForces()
      {  UTIL_THROW("Unimplemented method"); }

      /**
      * Calculate the total external energy for the associated System.
      */
      virtual double energy() const = 0;

   };

} 
#endif
#endif
