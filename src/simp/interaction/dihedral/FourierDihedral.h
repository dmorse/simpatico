#ifndef SIMP_FOURIER_DIHEDRAL_H
#define SIMP_FOURIER_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <util/containers/DArray.h>        // member
#include <util/space/Vector.h>             // inline function
#include <simp/interaction/dihedral/Torsion.h>        // inline function
#include <simp/interaction/dihedral/TorsionForce.h>   // inline function
#include <util/global.h>

namespace Simp
{

   using namespace Util;

   /**
   * A truncated Fourier series dihedral potential.
   *
   * This class defines a dihedral potential:
   * \f[
   *    V(\phi) = K_{0} + \sum_{m=1}^{4} K_{m} \cos(m \phi)]
   * \f]
   * where phi is the dihedral potential defined in \ref Simp_Interaction_Dihedral_Module.
   *
   * \sa \ref simp_interaction_dihedral_FourierDihedral_page "parameter file format"
   *   
   * \ingroup Simp_Interaction_Dihedral_Module
   */
   class FourierDihedral : public ParamComposite 
   {
   
   public:

      /**
      * Default constructor.
      */
      FourierDihedral();
   
      /**
      * Copy constructor.
      */
      FourierDihedral(const FourierDihedral& other);
   
      // Default C++ destructor.
 
      /**
      * Assignment.
      */
      FourierDihedral& operator = (const FourierDihedral& other);
  
      /**
      * Set the number of dihedral types.
      *
      * \param nDihedralType number of dihedral types
      */
      void setNDihedralType(int nDihedralType);
       
      /**
      * Read dihedral interaction parameters from input stream.
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
      * Returns potential energy for one dihedral group.
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
      * \param b1     bond vector from atom 1 to 2 (input)
      * \param b2     bond vector from atom 2 to 3 (input)
      * \param b3     bond vector from atom 3 to 4 (input)
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
      * \param name  parameter name (k0, k1, k2, k3, or k4)
      * \param type  dihedral type index
      * \return parameter value
      */
      double get(std::string name, int type) const;

      /**
      * Return name string "FourierDihedral".
      */
      std::string className() const;
 
   private:
  
      /// List of coefficients for one type of dihedral 
      class CoeffList;

      /// Array of coefficients, one struct per dihedral type
      DArray<FourierDihedral::CoeffList> coeffs_;

      /// Number of dihedral types
      int nDihedralType_;        

   // friends:

      friend std::istream& 
      operator >> (std::istream& in, FourierDihedral::CoeffList& param);

      friend std::ostream& 
      operator << (std::ostream& in, const FourierDihedral::CoeffList& param);

   };

   /**
   * List of coefficients for one type of dihedral.
   */
   class FourierDihedral::CoeffList 
   {

   public:

      // Coefficients in expansion V(phi) = \sum_m k_m cos(m phi)
      double k0; 
      double k1; 
      double k2;
      double k3; 
      double k4;

      // Coefficients in expansion V(phi) = \sum_m a_m cos^m(phi)
      double a0;
      double a1; 
      double a2; 
      double a3; 
      double a4; 

      // Coefficients in expansion of derivative dV/d(cos(phi))
      double g2; 
      double g3; 
      double g4; 

      // Default constructor.
      CoeffList() 
       : k0(0.0),
         k1(0.0),
         k2(0.0),
         k3(0.0),
         k4(0.0),
         a0(0.0),
         a1(0.0),
         a2(0.0),
         a3(0.0),
         a4(0.0),
         g2(0.0),
         g3(0.0),
         g4(0.0)
      {}

      void init()
      {
         a0 = k0 - k2 + k4;
         a1 = k1 - 3.0*k3;
         a2 = 2.0*k2 - 8.0*k4;
         a3 = 4.0*k3;
         a4 = 8.0*k4;

         g2 = 2.0*a2;
         g3 = 3.0*a3;
         g4 = 4.0*a4;
      }

      template <class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
         ar & k0; 
         ar & k1; 
         ar & k2;
         ar & k3; 
         ar & k4;
         ar & a0;
         ar & a1; 
         ar & a2; 
         ar & a3; 
         ar & a4; 
         ar & g2; 
         ar & g3; 
         ar & g4; 
      }

   };

   // Friend function declarations
 
   std::istream& 
   operator >> (std::istream& in, FourierDihedral::CoeffList& param);

   std::ostream& 
   operator << (std::ostream& in, const FourierDihedral::CoeffList& param);


   // Inline member function definitions
   
   /* 
   * Return dihedral energy.
   */
   inline
   double FourierDihedral::energy(const Vector& b1, const Vector& b2,
          const Vector& b3, int type) const
   {
      Torsion torsion;
      bool status;
      status = torsion.computeAngle(b1, b2, b3); // computes cosPhi
      if (!status) { 
         double c = torsion.cosPhi;
         const CoeffList* p = &coeffs_[type];
         return (p->a0 + c*(p->a1 + c*(p->a2 + c*(p->a3 + c*p->a4))));
      } else {
         return 0.0;
         std::cout << "!";
      }
   }

   /* 
   * Compute and return vectors:
   *
   *    f1 = d energy / d(b1)
   *    f2 = d energy / d(b2)
   *    f3 = d energy / d(b3)
   *
   * for use in MD and stress calculation.
   */
   inline 
   void FourierDihedral::force(const Vector& b1, const Vector& b2,
                const Vector& b3, Vector& f1, Vector& f2, Vector& f3, 
                int type) const
   {
      TorsionForce torsion;
      bool status;
      status = torsion.computeDerivatives(b1, b2, b3);
      if (!status) {
         double c = torsion.cosPhi;
         const CoeffList* p = &coeffs_[type];
         double dEdC = p->a1 + c*(p->g2 + c*(p->g3 + c*p->g4));
         f1.multiply(torsion.d1, dEdC);
         f2.multiply(torsion.d2, dEdC);
         f3.multiply(torsion.d3, dEdC);
      } else {
         f1.zero();
         f2.zero();
         f3.zero();
      }
   }

} 
#endif
