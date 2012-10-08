#ifndef INTER_TANH_COSINE_EXTERNAL_H
#define INTER_TANH_COSINE_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>
#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   * A clipped cosine potential that induces lamellar ordering
   * along the direction specified by perpDirection_.
   *    perpDirection_ = 0: x direction
   *                   = 1: y direction
   *                   = 2: z direction
   *
   *                                                  /                   /                   z   \ \
   *  u = prefactor[atomType] externalParameter tanh | clipParameter cos | 2  pi periodicity ---   | |
   *                                                  \                   \                   Lz  / / 
   *
   * Prefactor (which depends on the atomType), externalParameter, interfaceWidth (relative to the box length 
   * along the direction perpendicular to lamellae) and periodicity are given as inputs in the parameter file. 
   * ClipParameter is the inverse of 2*pi*periodicity*interfaceWidth. 
   *
   * \ingroup Inter_External_Module
   */
   class TanhCosineExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      TanhCosineExternal();

      /**
      * Copy constructor.
      */
      TanhCosineExternal(const TanhCosineExternal& other);

      /**
      * Assignment.
      */
      TanhCosineExternal& operator = (const TanhCosineExternal& other);

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Sets external parameter
      *
      * \param externalParameter external parameter of system
      */
      void setExternalParameter(double externalParameter);

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate length along perpendicular direction).
      */
      void setBoundary(Boundary &boundary);

      /**
      * Read potential parameters, and initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
      * \pre Boundary must have been set, by calling setBoundary().
      *
      * \param in input stream 
      */
      void readParameters(std::istream &in);

      /**
      * Returns external parameter
      *
      * \return external parameter
      */
      double externalParameter() const;

      /**
      * Returns external potential energy of a particle of type i.
      *
      * \param d  component of position along the perpendicular direction
      * \param i  type of particle (prefactor depends on atomtype)
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
      * Returns magnitude of the external force.
      *
      * \param d    component of position along the perpendicular direction
      * \param type atom type id (not used)
      * \return    force scalar
      */
      double forceScalar(double d, int type) const;
 
      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      void getForce(const Vector& position, int type, Vector& force) const;
 
      /**
      * Return name string "TanhCosineExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 2;

      /// Index representing the direction perpendicular to the lamellae.
      int perpDirection_;

      /// Interfcial width in lamellar phase.
      double width_;

      /// Prefactor array ofsize nAtomType.
      DArray<double> prefactor_;   

      /// External parameter.
      double externalParameter_;

      /// Number of periods in a cell
      int periodicity_;

      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 
 
    
   inline double TanhCosineExternal::energy(double d, int type) const
   {
      double perpLength_, q_, clipParameter_, arg, clipcos;
      Vector lengths_;
      lengths_ = boundaryPtr_->lengths();
      perpLength_ = lengths_[perpDirection_];
      q_ = (2.0*M_PI*periodicity_)/perpLength_;
      clipParameter_   = 1.0/(q_*width_*perpLength_);

      arg = q_*d;
      clipcos = clipParameter_*cos(arg);
      
      return prefactor_[type]*externalParameter_*tanh(clipcos);
   }

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline 
   double TanhCosineExternal::energy(const Vector& position, int type) const
   {
      double d, totalEnergy;
      totalEnergy = 0.0;
      d = position[perpDirection_];
      totalEnergy += energy(d, type);
      
      return totalEnergy;
   }

   /* 
   * Calculate force for a particle as a function of distance to boundary.
   */
   inline double TanhCosineExternal::forceScalar(double d, int type) const
   {
      double perpLength_, q_, clipParameter_, arg, clipcos, tanH, sechSq;
      Vector lengths_;
      lengths_ = boundaryPtr_->lengths();
      perpLength_ = lengths_[perpDirection_];
      q_ = (2.0*M_PI*periodicity_)/perpLength_;
      clipParameter_   = 1.0/(q_*width_*perpLength_);
      arg = q_*d;
      clipcos = clipParameter_*cos(arg);
      tanH = tanh(clipcos);
      sechSq = (1.0 - tanH*tanH);
      return prefactor_[type]*externalParameter_*sechSq*clipParameter_*sin(arg)*q_;
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void TanhCosineExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      double d = position[perpDirection_];
      double scalarf;
      force.zero();
      scalarf = forceScalar(d, type);
      force[perpDirection_] = scalarf;
   }
 
}
#endif
