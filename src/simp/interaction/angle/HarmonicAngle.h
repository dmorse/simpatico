#ifndef SIMP_HARMONIC_ANGLE_H
#define SIMP_HARMONIC_ANGLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <simp/interaction/angle/BendForce.h>      // used in inline function
#include <util/global.h>
#include <cmath>

namespace Util {class Random;}

namespace Simp
{

   using namespace Util;

   /**
   * A angle potential that is harmonic in the angle.
   *
   * This class implements an angle potential: 
   *
   *    V(theta) = 0.5  kappa ( theta - theta0 )^2.
   *
   * \sa \ref simp_interaction_angle_HarmonicAngle_page "parameter file format"
   *
   * \ingroup Simp_Interaction_Angle_Module
   */
   class HarmonicAngle : public ParamComposite
   {
   
   public:

      /**
      * Default constructor.
      */
      HarmonicAngle();
   
      /**
      * Copy constructor.
      */
      HarmonicAngle(const HarmonicAngle& other);
   
      /**
      * Assignment.
      */
      HarmonicAngle& operator = (const HarmonicAngle& other);
  
      // Default C++ destructor.
 
      /**
      * Set the number of angle types.
      *
      * \param nAngleType number of angle types
      */
      void setNAngleType(int nAngleType);
       
      /**
      * Read angle interaction parameters from input stream.
      *
      * Format:
      *     kappa   CArray<double>
      *     theta0  CArray<double>
      *
      * \pre nAngleType must be set, by calling setNAngleType().
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
      * Modify a parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  angle type index 
      * \param value new value of parameter
      */
      void set(std::string name, int type, double value);

      /**
      * Returns potential energy for one angle.
      *
      * \param cosTheta  cosine of the bend angle.
      * \param type      type of bend angle.
      */
      double energy(double cosTheta, int type) const;
 
      /**
      * Compute angle forces.
      *
      * Computes forces F1 and F2 along the two bonds in the angle.
      *
      * \param b1     bond vector from atom 1 to 2.
      * \param b2     bond vector from atom 2 to 3.
      * \param F1     return force along b1 direction.
      * \param F2     return force along b2 direction.
      * \param type   type of angle.
      */
      void force(const Vector& b1, const Vector& b2,
                       Vector& F1, Vector& F2, int type) const;

      /**
      * Return bond angle chosen from equilibrium distribution.
      *
      * This function returns a bond angle chosen from the Boltzmann
      * distribution of angle for bonds of random orientation. The
      * distribution P(theta) of values of the angle theta is proportional 
      * to sin(theta)*exp[-beta*phi(theta) ], where phi(theta) 
      * is the bond energy. 
      *
      * \param random pointer to random number generator object.
      * \param beta   inverse temperature
      * \param type   bond type
      * \return random bond length chosen from equilibrium distribution.
      */
      double randomAngle(Random *random, double beta, int type) const;

      /**
      * Return bond angle cosine chosen from equilibrium distribution.
      *
      * This function returns a bond angle chosen from the Boltzmann
      * distribution of angle for bonds of random orientation. The
      * distribution P(theta) of values of the angle theta is proportional 
      * to sin(theta)*exp[-beta*phi(theta) ], where phi(theta) 
      * is the bond energy. 
      *
      * \param random pointer to random number generator object.
      * \param beta   inverse temperature
      * \param type   bond type
      * \return random bond length chosen from equilibrium distribution.
      */
      double randomCosineAngle(Random *random, double beta, int type) const;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param  name  parameter name
      * \param  type  angle type index
      * \return parameter value
      */
      double get(std::string name, int type) const;

      /**
      * Return name string "HarmonicAngle" for this evaluator class.
      */
      std::string className() const;
 
   private:
   
      /// Maximum possible number of angle types (two bonds: 11, 12, 21, 22).
      static const int MaxNAngleType = 4;
   
      double  kappa_[MaxNAngleType];    ///< spring constant
      double  theta0_[MaxNAngleType];   ///< preferred angle
      int  nAngleType_;                 ///< number of angle types

   };
   
   // Inline method definitions
   
   /*
   * Return angle energy.
   */
   inline double HarmonicAngle::energy(double cosTheta, int type) const
   {
      double dTheta = std::acos(cosTheta) - theta0_[type];
      return 0.5*kappa_[type]*dTheta*dTheta;
   }

   /* 
   * Return:
   *    F1 = d energy / d(b1)
   *    F2 = d energy / d(b2)
   * for use in MD and stress calculation.
   */
   inline
   void HarmonicAngle::force(const Vector& b1, const Vector& b2,
                             Vector& F1, Vector& F2, int type) const
   {
      BendForce bend;
      bend.computeDerivatives(b1, b2);
      double s = bend.sinTheta();
      if (s > 1.0E-10) {
         double factor = kappa_[type]*(theta0_[type] - bend.theta())/s;
         F1.multiply(bend.d1, factor);
         F2.multiply(bend.d2, factor);
      } else {
         F1.zero();
         F2.zero();
      }
   }

}
 
#endif
