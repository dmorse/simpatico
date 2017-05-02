#ifndef SIMP_COSINE_DIHEDRAL_H
#define SIMP_COSINE_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/space/Vector.h>
#include <util/param/ParamComposite.h>  // base class

#include "Torsion.h"                    // used in in-line function
#include "TorsionForce.h"               // used in in-line function

//#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A dihedral potential proportional to cos(phi).
   *
   * This class defines a dihedral potential:
   * \f[
   *    V(\phi) = kappa [1 + cos(phi)]
   * \f]
   * where phi is the dihedral potential, as defined in \ref Simp_Interaction_Dihedral_Module.
   *
   * \sa \ref simp_interaction_dihedral_CosineDihedral_page "parameter file format"
   *
   * \ingroup Simp_Interaction_Dihedral_Module
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
      *     0   2    3
      *     o   o----o
      *      \ /
      *       o 
      *       1 
      *
      * \param b1     bond vector from atom 0 to 1
      * \param b2     bond vector from atom 1 to 2
      * \param b3     bond vector from atom 2 to 3
      * \param type   type id for dihedral group
      */
      double energy(const Vector& b1, const Vector& b2, const Vector& b3,
          int type) const;
 
      /**
      * Returns derivatives of energy with respect to bond vectors forming the
      * dihedral group.
      *
      * \param b1     bond vector from atom 1 to 2.
      * \param b2     bond vector from atom 2 to 3.
      * \param b3     bond vector from atom 3 to 4.
      * \param f1     derivative of energy w/respect to b1 (output)
      * \param f2     derivative of energy w/respect to b2 (output)
      * \param f3     derivative of energy w/respect to b3 (output)
      * \param type   type id for dihedral group
      */
      void force(const Vector& b1, const Vector& b2, const Vector& b3,
                 Vector& f1, Vector& f2, Vector& f3, int type) const;

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
   double CosineDihedral::energy(const Vector& b1, const Vector& b2,
          const Vector& b3, int type) const
   {
      Torsion torsion;
      bool status;
      status = torsion.computeAngle(b1, b2, b3); // computes cosPhi
      if (!status) {
         return (kappa_[type] * (1.0 + torsion.cosPhi));
      } else {
         return 0.0;
      }
   }

   /* 
   * Return:
   *    f1 = d energy / d(b1)
   *    f2 = d energy / d(b2)
   *    f3 = d energy / d(b3)
   * for use in MD and stress calculation.
   */
   inline
   void CosineDihedral::force(const Vector& b1, const Vector& b2,
        const Vector& b3, Vector& f1, Vector& f2, Vector& f3, int type) const
   {
      TorsionForce torsion;
      bool status; // Error code, 0 is normal, 1 is error
      status = torsion.computeDerivatives(b1, b2, b3);

      if (!status) {
         double dEdCosPhi = kappa_[type];
         f1.multiply(torsion.d1, dEdCosPhi);
         f2.multiply(torsion.d2, dEdCosPhi);
         f3.multiply(torsion.d3, dEdCosPhi);
      } else {
         f1.zero();
         f2.zero();
         f3.zero();
      }
   }

} 
#endif
