#ifndef INTER_MULTI_HARMONIC_DIHEDRAL_H
#define INTER_MULTI_HARMONIC_DIHEDRAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <util/space/Vector.h>             // inline function
#include <inter/dihedral/Torsion.h>        // inline function
#include <inter/dihedral/TorsionForce.h>   // inline function
#include <util/global.h>

namespace Inter
{

   using namespace Util;

   /**
   * A dihedral potential proportional to cos(phi).
   *
   * This class defines a dihedral potential:
   * \f[
   *    V(\phi) = K_{0} + \sum_{m=1}^{4} K_{m} \cos(m \phi)]
   * \f]
   * where phi is the dihedral potential, as defined in \ref Inter_Dihedral_Module.
   *
   * \ingroup Inter_Dihedral_Module
   */
   class MultiHarmonicDihedral : public ParamComposite 
   {
   
   public:

      /**
      * Default constructor.
      */
      MultiHarmonicDihedral();
   
      /**
      * Copy constructor.
      */
      MultiHarmonicDihedral(const MultiHarmonicDihedral& other);
   
      // Default C++ destructor.
 
      /**
      * Assignment.
      */
      MultiHarmonicDihedral& operator = (const MultiHarmonicDihedral& other);
  
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
      * Return name string "MultiHarmonicDihedral" for this evaluator class.
      */
      std::string className() const;
 
   private:
   
      struct Parameter 
      {

         // Coefficients in expansion V(phi) = \sum_m K_m cos(m phi)
         double K0; 
         double K1; 
         double K2;
         double K3; 
         double K4;

         // Coefficients in expansion V(phi) = \sum_m A_m cos^m(phi)
         double A0; 
         double A1; 
         double A2; 
         double A3; 
         double A4; 

         void init()
         {
            A0 = K0 - K2 + K4;
            A1 = K1 - 3.0*K3;
            A2 = 2.0*K2 - 8.0*K4;
            A3 = 4.0*K3;
            A4 = 8.0*K4;
         }

      }

      /// Array of parameters, one struct per dihedral type
      DArray<Parameters> parameters;

      /// Number of dihedral types
      int nDihedralType_;        

   };
 
   // Inline method definitions
   
   /* 
   * Return dihedral energy.
   */
   inline
   double MultiHarmonicDihedral::energy(const Vector& b1, const Vector& b2,
          const Vector& b3, int type) const
   {
      Torsion torsion;
      torsion.computeAngle(b1, b2, b3); // computes cosPhi
      double c = torsion.cosPhi;

      Parameter* p = &parameters[type];
      return p->A0 + c*(p->A1 + c*(p->A2 + c*(p->A3 + c*p->A4)));
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
   void MultiHarmonicDihedral::force(const Vector& b1, const Vector& b2,
        const Vector& b3, Vector& f1, Vector& f2, Vector& f3, int type) const
   {
      TorsionForce torsion;
      torsion.computeDerivatives(b1, b2, b3);
      double c = torsion.cosPhi;

      Parameter* p = parameters[type];
      double dEdC = p.A1 + c*(2.0*p.A2 + c*(3.0*p.A3 + 4.0*c*p.A4));
      f1.multiply(torsion.d1, dEdC);
      f2.multiply(torsion.d2, dEdC);
      f3.multiply(torsion.d3, dEdC);
   }

} 
#endif
