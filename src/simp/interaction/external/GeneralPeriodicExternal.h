#ifndef SIMP_GENERAL_PERIODIC_EXTERNAL_H
#define SIMP_GENERAL_PERIODIC_EXTERNAL_H

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
   *                                                 /              /                    /     /      w_i.(x-shifts_[0])    w_i.(y-shifts_[1])    w_i.(z-shifts_[2])            \  \  \ \
   * u = prefactor[atomType] externalParameter tanh | clipParameter|  Sum Prefactor_[i] | cos | 2 pi ------------------ + ------------------ + ------------------ + phase_i  |  |  | |
   *                                                 \              \  i                 \     \             Lx               Ly                       Lz                    /  /  / /
   *
   * Prefactor (which depends on the atomType), externalParameter, waveIntVectors, interfaceWidth and periodicity
   * are given as inputs in the parameter file. 
   * ClipParameter is the inverse of 2*pi*periodicity*interfaceWidth. 
   *
   * \ingroup Simp_Interaction_External_Module
   */
   class GeneralPeriodicExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      GeneralPeriodicExternal();

      /**
      * Copy constructor.
      */
      GeneralPeriodicExternal(const GeneralPeriodicExternal& other);

      /**
      * Assignment.
      */
      GeneralPeriodicExternal& operator = (const GeneralPeriodicExternal& other);

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
      * Return name string "GeneralPeriodicExternal".
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
      DArray<double> shifts_;

      /// Number of unit cells in box
      int periodicity_;

      /// Interface width
      double interfaceWidth_;

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
   inline double GeneralPeriodicExternal::energy(const Vector& position, int type) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*periodicity_*interfaceWidth_);
      
      Vector r;
      Vector a;
      r.zero();
      a.zero();
      double cosine = 0.0;
      for (int i = 0; i < nWaveVectors_; ++i) {
         r = position;
         for (int j = 0; j < Dimension; j++) {
            a[j] = shifts_[i*Dimension + j];
         }
         r -= a;
         Vector q;
         q[0] = 2.0*M_PI*periodicity_*waveVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*periodicity_*waveVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*periodicity_*waveVectors_[i][2]/cellLengths[2];
         double arg = q.dot(r)+phases_[i];
         cosine += prefactor_[i*nAtomType_ + type]*cos(arg);
      }
      cosine *= clipParameter;
      return externalParameter_*tanh(cosine);
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void GeneralPeriodicExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*periodicity_*interfaceWidth_);
 
      Vector r;
      Vector a;
      r.zero();
      a.zero();
      double cosine = 0.0;
      Vector deriv;
      deriv.zero();
      for (int i = 0; i < nWaveVectors_; ++i) {
         r = position;
         for (int j = 0; j < Dimension; j++) {
            a[j] = shifts_[i*Dimension + j];
         }
         r -= a;
         Vector q;
         q[0] = 2.0*M_PI*periodicity_*waveVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*periodicity_*waveVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*periodicity_*waveVectors_[i][2]/cellLengths[2];
         double arg = q.dot(r)+phases_[i];
         cosine += prefactor_[i*nAtomType_ + type]*cos(arg);
         double sine = -1.0*sin(arg);
         q *= prefactor_[i*nAtomType_ + type]*sine;
         deriv += q;
      }
      cosine *= clipParameter;
      deriv *= clipParameter;
      double tanH = tanh(cosine);
      double sechSq = (1.0 - tanH*tanH);
      double f = externalParameter_*sechSq;
      deriv *= -1.0*f;
      force = deriv;
   }
 
}
#endif
