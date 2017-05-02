#ifndef SIMP_BOX_EXTERNAL_H
#define SIMP_BOX_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>
#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A repulsive 9-3 potential confining particles within a box.
   *
   *      4 epsilon pi rho   / sigma^9        sigma^6       25     \
   *  u = ----------------- | --------- - 7.5 -------- + ---------- |
   *       45 sigma^(-3)     \   r^9            r^6       sqrt(10) /
   *
   * The coefficients are obtained by integrating the repulsive LJPair
   * potential in a half space with default density rho = 0.7. Member variables
   * include strength of the potential epsilon, the range of interaction sigma,
   * and the cutoff length 0.4^(1/6) = 0.85837422. Member functions evaluate
   * energy and force for individual particles. The internal cutoff 0.506079 =
   * 0.129615^(1/3) is set to yield a cutoff potential of 80 for the default
   * density 0.7.
   *
   * The potential is independent on the type of atoms.
   * 
   * \ingroup Simp_Interaction_External_Module
   */
   class BoxExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      BoxExternal();

      /**
      * Copy constructor.
      */
      BoxExternal(const BoxExternal& other);

      /**
      * Assignment.
      */
      BoxExternal& operator = (const BoxExternal& other);

      /// \name Mutators
      //@{ 

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate distances).
      */
      void setBoundary(Boundary &boundary);

      /**
      * Sets external parameter
      *
      * \param externalParameter external parameter of system
      */
      void setExternalParameter(double externalParameter);

      /**
      * Read potential parameters, and initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
      * \pre Boundary must have been set, by calling setBoundary().
      *
      * \param in  input stream 
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

      /*
      * Set a potential energy parameter, identified by a string.
      */
      void set(std::string name, double value);

      /*
      * Get a parameter value, identified by a string.
      */
      double get(std::string name) const;

      //@}

      /// \name Accessors
      //@{ 

      /**
      * Return name string "LJPair" for this evaluator class.
      */
      std::string className() const;

      /**
      * Returns external parameter
      *
      * \return external parameter
      */
      double externalParameter() const;
 
      /**
      * Returns external potential energy of a particle of type i.
      *
      * \param d  distance to the nearest boundary
      * \param i  type of particle (not used)
      * \return   external potential energy
      */
      double energy(double d, int i) const;
 
      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return     external potential energy
      */
      double energy(const Vector& position, int i) const;

      /**
      * Returns magnitude of the repulsive force.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector.
      *
      * \param d   distance from the nearest wall.
      * \param i   atom type id (not used)
      * \return    force scalar
      */
      double forceScalar(double d, int i) const;
 
      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      void getForce(const Vector& position, int type, Vector& force) const;
 
      //@}

   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 2;

      /// Strength of the confining potential.
      double epsilon_;

      /// Range of the repulsion.
      double sigma_;

      /// Range of the repulsion squared.
      double sigmaCb_;
 
      /// Cutoff distance.
      double cutoff_;

      /// Prefactor.
      double coeff_;
  
      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 
 
   /* 
   * Calculate repulsive potential energy for particle of type i.
   *
   * 7.905694 = 25/sqrt(10), which is the negative of the minimum of the
   * unshifted 9-3 potential.
   */
   inline double BoxExternal::energy(double d, int type) const
   {
      if (d < cutoff_) {
         double r3i;
         r3i = d * d * d;
         if (r3i < 0.129615*sigmaCb_) r3i = 0.129615 * sigmaCb_;
         r3i = sigmaCb_ / r3i;
         return coeff_ * (r3i*(r3i*r3i - 7.5) + 7.905694);
      } else {
         return 0.0;
      }
   }
   
   /* 
   * Calculate external potential energy for a single atom.
   */
   inline double BoxExternal::energy(const Vector& position, int type) const
   {
      double d, totalEnergy;
      Vector halfL = boundaryPtr_->lengths();
      halfL *= 0.5;
      totalEnergy = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         d = position[i];
         if (d > halfL[i]) d = halfL[i] + halfL[i] - d;
            totalEnergy += energy(d, type);
      }
      return totalEnergy;
   }
   
   /* 
   * Calculate force for a particle as a function of distance to boundary.
   */
   inline double BoxExternal::forceScalar(double d, int i) const
   {
      if (d < cutoff_) {
         double r3i;
         r3i = d * d * d;
         if (r3i < 0.129615 * sigmaCb_) r3i = 0.129615 * sigmaCb_;
         r3i = sigmaCb_ / r3i;
         return coeff_ * (9.0*r3i*r3i - 22.5) * r3i / d;
      } else {
         return 0.0;
      }
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline void 
   BoxExternal::getForce(const Vector& position, int type, Vector& force) const
   {
      double d = position[2];
      double scalarf;
      Vector halfL = boundaryPtr_->lengths();
      halfL *= 0.5;
      for (int i = 0; i < Dimension; ++i) {
         d = position[i];
         if (d > halfL[i]) {
            d = halfL[i] + halfL[i] - d;
            scalarf = forceScalar(d, type);
            force[i] = -scalarf;
         } else {
            scalarf = forceScalar(d, type);
            force[i] = scalarf;
         }
      }
   }
 
}
#endif
