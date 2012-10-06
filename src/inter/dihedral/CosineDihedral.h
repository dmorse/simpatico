#ifndef INTER_COSINE_DIHEDRAL_H
#define INTER_COSINE_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/space/Vector.h>
#include <util/param/ParamComposite.h>  // base class

#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   * A four-body dihedral potential.
   *
   * Simple dihedral potential: kappa [1 - cos(theta)], where theta
   * is the dihedral angle.
   *
   * \ingroup Inter_Dihedral_Module
   */
   class CosineDihedral : public ParamComposite 
   {
   
   public:

      /**
      * Default constructor.
      */
      CosineDihedral();
   
      /**
      * Copy constructor.
      */
      CosineDihedral(const CosineDihedral& other);
   
      /**
      * Assignment.
      */
      CosineDihedral& operator = (const CosineDihedral& other);
  
      // Default C++ destructor.
 
      /**
      * Set the number of dihedral types.
      *
      * \param nDihedralType number of dihedral types
      */
      void setNDihedralType(int nDihedralType);
       
      /**
      * Read dihedral interaction parameters from input stream.
      *
      * Format:
      *     kappa  CArray<double>
      *
      * \pre nDihedralType must be set, by calling setNDihedralType().
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
      * \param type  dihedral type index 
      * \param value new value of parameter
      */
      void set(std::string name, int type, double value);

      /**
      * Returns potential energy for one dihedral.
      *
      *     1   3    4
      *     o   o----o
      *      \ /
      *       o 
      *       2 
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param type   type of dihedral.
      */
      double energy(const Vector& R1, const Vector& R2, const Vector& R3,
          int type) const;
 
      /**
      * Returns derivatives of energy with respect to bond vectors forming the
      * dihedral group.
      *
      * \param R1     bond vector from atom 1 to 2.
      * \param R2     bond vector from atom 2 to 3.
      * \param R3     bond vector from atom 3 to 4.
      * \param F1     return force along R1 direction.
      * \param F2     return force along R2 direction.
      * \param F3     return force along R2 direction.
      * \param type   type of dihedral.
      */
      void force(const Vector& R1, const Vector& R2, const Vector& R3,
                 Vector& F1, Vector& F2, Vector& F3, int type) const;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name  parameter name
      * \param type  dihedral type index
      * \return parameter value
      */
      double get(std::string name, int type) const;

      /**
      * Return name string "CosineDihedral" for this evaluator class.
      */
      std::string className() const;
 
   private:
   
      /// Maximum possible number of dihedral types (assuming 2 types of bonds).
      static const int MaxNDihedralType = 8;
   
      double kappa_[MaxNDihedralType];  ///< spring constant
      int    nDihedralType_;            ///< number of dihedral types

   };
   
   // Inline method definitions
   
   /* 
   * Return dihedral energy.
   */
   inline
   double CosineDihedral::energy(const Vector& R1, const Vector& R2,
          const Vector& R3, int type) const
   {
      Vector u1, u2;

      u1.cross(R1, R2);
      u1 /= u1.abs();

      u2.cross(R2, R3);
      u2 /= u2.abs();

      return ( kappa_[type] * (1.0 + u1.dot(u2)) );
   }

   /* 
   * Return:
   *    F1 = d energy / d(R1)
   *    F2 = d energy / d(R2)
   *    F3 = d energy / d(R3)
   * for use in MD and stress calculation.
   */
   inline
   void CosineDihedral::force(const Vector& R1, const Vector& R2,
        const Vector& R3, Vector& F1, Vector& F2, Vector& F3, int type) const
   {
      Vector u1, u2, tmp1, tmp2;
      double r1, r2, cosPhi;

      u1.cross(R1, R2);
      r1 = u1.abs();
      u1 /= r1;

      u2.cross(R2, R3);
      r2 = u2.abs();
      u2 /= r2;

      cosPhi = u1.dot(u2);

      tmp1.multiply(u1, -cosPhi);
      tmp1 += u2;
      tmp1 *= kappa_[type] / r1;

      tmp2.multiply(u2, -cosPhi);
      tmp2 += u1;
      tmp2 *= kappa_[type] / r2;

      F1.cross(R2, tmp1);

      F2.cross(tmp1, R1);
      tmp1.cross(R3, tmp2);
      F2 += tmp1;

      F3.cross(tmp2, R2);
   }

} 
#endif
