#ifdef  MCMD_TETHER
#ifndef HARMONIC_TETHER_H
#define HARMONIC_TETHER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>

#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A harmonic tether potential.
   *
   * This class implements a harmonic potential with a zero rest length.
   *
   * \ingroup Simp_Interaction_Tether_Module
   */
   class HarmonicTether : public ParamComposite 
   {

   public:
 
      /**
      * Default constructor.
      */
      HarmonicTether();

      /**
      * Copy constructor.
      */
      HarmonicTether(const HarmonicTether& other);

      /**
      * Assignment.
      */
      HarmonicTether& operator = (const HarmonicTether& other);

      // Default constructor is fine.

      /**
      * Read tether interaction parameters from file.
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

      /**
      * Returns potential energy for one tether.
      *
      * \param rSq  square of distance between tethered particles.
      * \param type type of tether.
      */
      double energy(double rSq, int type) const;
  
      /**
      * Returns force/distance for one tether, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between atom and anchor.
      * \param type type of tether.
      * \return     scalar repulsive force divided by distance.
      */
      double forceOverR(double rSq, int type) const;
  
   private:
   
      /// Maximum possible number of tether types
      static const int MaxNTetherType = 2;
   
      /// Spring constants.
      double kappa_[MaxNTetherType];  

      /// Number of tether types
      int    nTetherType_;            

   };
   
   // Inline method definitions
   
   /* 
   * Return interaction energy for one tether
   */
   inline double HarmonicTether::energy(double rSq, int type) const
   {  return 0.5*kappa_[type]*rSq; }

   /* 
   * Return force / distance for one tether, for use in MD
   */
   inline double HarmonicTether::forceOverR(double rSq, int type) const
   {  return -kappa_[type]; }
  
} 
#endif
#endif
