#ifndef SIMP_HARMONIC_L0_BOND_H
#define SIMP_HARMONIC_L0_BOND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <util/param/ParamComposite.h>  // base class
#include <util/random/Random.h>           // Util namespace

#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A harmonic bond potential with zero rest length.
   *
   * \sa \ref simp_interaction_bond_HarmonicL0Bond_page "Parameter file format"
   * \sa \ref simp_interaction_bond_interface_page
   *
   * \ingroup Simp_Interaction_Bond_Module
   */
   class HarmonicL0Bond : public ParamComposite 
   {
   
   public:

      /**
      * Constructor.
      */
      HarmonicL0Bond();
   
      /**
      * Copy constructor.
      */
      HarmonicL0Bond(const HarmonicL0Bond& other);
   
      // Default destructor.
 
      /**
      * Assignment.
      */
      HarmonicL0Bond& operator = (const HarmonicL0Bond& other);
  
      /**
      * Set the number of bond types.
      *
      * \param nBondType number of bond types
      */
      void setNBondType(int nBondType);

      /**
      * Read bond interaction parameters from input stream.
      *
      * Format:
      *     kappa  CArray<double>
      *
      * \pre nBondType must be set, by calling setNBondType().
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

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

      /**
      * Returns potential energy for one bond.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      */
      double energy(double rSq, int type) const;
   
      /**
      * Returns force/distance for one bond, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      * \return     scalar repulsive force divided by distance.
      */
      double forceOverR(double rSq, int type) const;
   
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
      * \return random bond length chosen from equilibrium distribution.
      */
      double randomBondLength(Random *random, double beta, int type) const;
   
      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param type   bond type index 
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value);

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter name
      * \param type  bond type index
      */
      double get(std::string name, int type) const;

      /**
      * Returns bond spring constant
      */
      double kappa(int type) const;

   private:
   
      /// Maximum possible number of bond types
      static const int MaxNBondType = 2;
   
      double kappa_[MaxNBondType];  ///< spring constant
      int    nBondType_;            ///< number of bond types

   };
   
   // Inline method definitions

   /*
    * Return bond spring constant
    */
   inline double HarmonicL0Bond::kappa(int type) const
   {
      assert(type >= 0);
      assert(type < nBondType_);
      return kappa_[type];
   }
   
   /* 
   * Return interaction energy for one bond
   */
   inline double HarmonicL0Bond::energy(double rSq, int type) const
   {  return 0.5*kappa_[type]*rSq; }

   /* 
   * Return force / distance for one bond, for use in MD
   */
   inline 
   double HarmonicL0Bond::forceOverR(double rSq, int type) const
   {  return -kappa_[type]; }
   
} 
#endif
