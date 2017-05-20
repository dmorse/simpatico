#ifndef SIMP_NUCLEATION_EXTERNAL_H
#define SIMP_NUCLEATION_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>
#include <cmath>

namespace Simp
{

   using namespace Util;

   /**
   * A clipped cosine potential that induces ordering
   * along directions specified by waveIntVectors, w_i.
   *
   *                                                 /              /          /     /      w_i.(x-shift_[0])    w_i.(y-shift_[1])    w_i.(z-shift_[2])            \  \  \ \
   * u = prefactor[atomType] externalParameter tanh | clipParameter| C_ +  Sum | cos | 2 pi ------------------ + ------------------ + ------------------ + phase_i  |  |  | |
   *                                                 \              \      i   \     \             Lx               Ly                       Lz                    /  /  / /
   *
   * Prefactor (which depends on the atomType), externalParameter, waveIntVectors, interfaceWidth and periodicity
   * are given as inputs in the parameter file. 
   * ClipParameter is the inverse of 2*pi*periodicity*interfaceWidth. 
   *
   * \ingroup Simp_Interaction_External_Module
   */
   class NucleationExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      NucleationExternal();

      /**
      * Copy constructor.
      */
      NucleationExternal(const NucleationExternal& other);

      /**
      * Assignment.
      */
      NucleationExternal& operator = (const NucleationExternal& other);

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

      /**
      * Returns external parameter
      *
      * \return external parameter
      */
      double externalParameter() const;

      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return     external potential energy
      */
      double energy(const Vector& position, int i) const;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      void getForce(const Vector& position, int type, Vector& force) const;
 
      /**
      * Return name string "NucleationExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 3;

      /// Prefactor array ofsize nAtomType.
      DArray<double> prefactor_;

      /// External parameter.
      double externalParameter_;

      /// Number of reciprocal lattice vectors
      int  nWaveVectors_;

      /// Array of Miller index IntVectors for the reciprocal lattice vectors.
      DArray<Vector>  waveVectors_;

      /// Phases for the different plane waves.
      DArray<double> phases_;

      /// Prefactor array ofsize nAtomType.
      Vector shift_;

      /// Prefactor array ofsize nAtomType.
      double C_;

      /// Number of unit cells in box
      int periodicity_;

      /// Interface width
      double interfaceWidth_;

      /// Prefactor array ofsize nAtomType.
      double nucleationClip_;

      /// Prefactor array ofsize nAtomType.
      double bias_;

      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline double NucleationExternal::energy(const Vector& position, int type) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*periodicity_*interfaceWidth_);
      double e;
      
      Vector r = position;
      r -= shift_;
      double cosine = 0.0;
      for (int i = 0; i < nWaveVectors_; ++i) {
         Vector q;
         q[0] = 2.0*M_PI*periodicity_*waveVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*periodicity_*waveVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*periodicity_*waveVectors_[i][2]/cellLengths[2];
         double arg = q.dot(r)+phases_[i];
         cosine += cos(arg);
      }

      cosine *= clipParameter;
      e = prefactor_[type]*externalParameter_*tanh(C_+cosine) * 
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))))+1)/2 *
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))))+1)/2 * 
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))))+1)/2 ;

      return e;
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void NucleationExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*periodicity_*interfaceWidth_);
      double ne;

      ne = 
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))))+1)/2 * 
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))))+1)/2 * 
      (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))))+1)/2 ;

      Vector r = position;
      r -= shift_;
      double cosine = 0.0;
      Vector deriv;
      deriv.zero();
      for (int i = 0; i < nWaveVectors_; ++i) {
         Vector q;
         q[0] = 2.0*M_PI*periodicity_*waveVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*periodicity_*waveVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*periodicity_*waveVectors_[i][2]/cellLengths[2];
         double arg = q.dot(r)+phases_[i];
         cosine += cos(arg);
         double sine = -1.0*sin(arg);
         q *= sine;
         deriv += q;
      }
      cosine *= clipParameter;

      deriv *= clipParameter;
      double tanH = tanh(C_+cosine);
      double sechSq = (1.0 - tanH*tanH);
      double f = prefactor_[type]*externalParameter_*sechSq;
      deriv *= -1.0*f;
      deriv *= ne;
      deriv[0] = deriv[0] + nucleationClip_*(M_PI/cellLengths[0]*sin(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_)))/
                            cosh(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))/
                            cosh(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))))+1)/2*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))))+1)/2;

      deriv[1] = deriv[1] + nucleationClip_*(M_PI/cellLengths[1]*sin(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_)))/
                            cosh(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))/
                            cosh(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))))+1)/2*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))))+1)/2;

      deriv[2] = deriv[2] + nucleationClip_*(M_PI/cellLengths[2]*sin(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_)))/
                            cosh(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))/
                            cosh(2.0*M_PI*position[2]/cellLengths[2]+acos(bias_))*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[0]/cellLengths[0]+acos(bias_))))+1)/2*
                            (tanh(nucleationClip_*(-bias_+cos(2.0*M_PI*position[1]/cellLengths[1]+acos(bias_))))+1)/2;
      force = deriv;
   }
 
}
#endif
