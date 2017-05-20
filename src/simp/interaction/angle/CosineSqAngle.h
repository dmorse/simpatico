#ifndef SIMP_COSINE_SQ_ANGLE_H
#define SIMP_COSINE_SQ_ANGLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/space/Vector.h>
#include <util/param/ParamComposite.h>  // base class

#include <cmath>

namespace Util {class Random;}

namespace Simp
{

   using namespace Util;

   /**
   * A three body angle potential, as a function of angle cosine.
   *
   * This class implements an angle potential: 
   *    kappa (cos[theta0] - cos[theta])^2.
   *
   * \sa \ref simp_interaction_angle_CosineSqAngle_page "parameter file format"
   *
   * \ingroup Simp_Interaction_Angle_Module
   */
   class CosineSqAngle : public ParamComposite 
   {
   
   public:

      /**
      * Default constructor.
      */
      CosineSqAngle();
   
      /**
      * Copy constructor.
      */
      CosineSqAngle(const CosineSqAngle& other);
   
      /**
      * Assignment.
      */
      CosineSqAngle& operator = (const CosineSqAngle& other);
  
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
      *     kappa  CArray<double>
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
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param type   type of angle.
      */
      void force(const Vector& R1, const Vector& R2,
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
      * Return name string "CosineSqAngle" for this evaluator class.
      */
      std::string className() const;
 
   private:
   
      /// Maximum possible number of angle types (two bonds: 11, 12, 21, 22).
      static const int MaxNAngleType = 4;
   
      double kappa_[MaxNAngleType];      ///< spring constant
      double theta0_[MaxNAngleType];     ///< spring constant
      double cosTheta0_[MaxNAngleType];  ///< spring constant
      int    nAngleType_;                ///< number of angle types

   };
   
   // Inline method definitions
   
   /* 
   * Return angle energy.
   */
   inline double CosineSqAngle::energy(double cosTheta, int type) const
   {
      return ( kappa_[type]
              *(cosTheta-cosTheta0_[type])
              *(cosTheta-cosTheta0_[type])*0.5 );
   }

   /* 
   * Return:
   *    F1 = d energy / d(R1)
   *    F2 = d energy / d(R2)
   * for use in MD and stress calculation.
   */
   inline
   void CosineSqAngle::force(const Vector& R1, const Vector& R2,
                   Vector& F1, Vector& F2, int type) const
   {
      Vector u1 = R1;
      Vector u2 = R2;
      double r1 = R1.abs();
      double r2 = R2.abs();
      double cosTheta;
     
      u1 /= r1;
      u2 /= r2;
      cosTheta = u1.dot(u2);

      F1.multiply(u1, cosTheta);
      F1 -= u2;
      F1 *= kappa_[type]/r1*(cosTheta0_[type]-cosTheta);

      F2.multiply(u2, cosTheta);
      F2 -= u1;
      F2 *= kappa_[type]/r2*(cosTheta0_[type]-cosTheta);
   }
  
}
 
#endif
