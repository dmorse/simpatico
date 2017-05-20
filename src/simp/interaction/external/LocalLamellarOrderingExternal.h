#ifndef SIMP_LOCAL_LAMELLAR_ORDERING_EXTERNAL_H
#define SIMP_LOCAL_LAMELLAR_ORDERING_EXTERNAL_H

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
   * A clipped cosine potential that induces lamellar ordering
   * along the direction specified by perpDirection_ in only a
   * fraction of box specified by the parameter fraction along
   * box length in parallelDirection_.
   *    perpDirection_ = 0: x direction
   *                   = 1: y direction
   *                   = 2: z direction
   *    parallelDirection_ = 0: x direction
   *                       = 1: y direction
   *                       = 2: z direction
   *
   *                                                  /                   /                   z   \ \
   *  u = prefactor[atomType] externalParameter tanh | clipParameter cos | 2  pi periodicity ---   | |
   *                                                  \                   \                   Lz  / / 
   *                                                  
   *                                   /          /      x - 0.5*(1-f)*Lx  \          /      0.5*(1+f)*Lx - x \  \
   *                           * 0.5* | 1 + tanh | 2 pi ------------------  | * tanh | 2 pi ------------------ |  |
   *                                   \          \            f*Lx        /          \             f*Lx      /  /
   *
   * Prefactor (which depends on the atomType), externalParameter, interfaceWidth (relative to the box length 
   * along the direction perpendicular to lamellae) and periodicity are given as inputs in the parameter file. 
   * ClipParameter is the inverse of 2*pi*periodicity*interfaceWidth. 
   *
   * \ingroup Simp_Interaction_External_Module
   */
   class LocalLamellarOrderingExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      LocalLamellarOrderingExternal();

      /**
      * Copy constructor.
      */
      LocalLamellarOrderingExternal(const LocalLamellarOrderingExternal& other);

      /**
      * Assignment.
      */
      LocalLamellarOrderingExternal& operator = (const LocalLamellarOrderingExternal& other);

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
      * \param position  atomic position Vector
      * \param i  atom type.
      * \return  external potential energy
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
      * Return name string "LocalLamellarOrderingExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 2;

      /// Index representing the direction perpendicular to the lamellae.
      int perpDirection_;

      /// Index representing the direction parallel to the lamellae.
      int parallelDirection_;

      /// Fraction of length in parallel direction for which lamellar ordering is imposed.
      double fraction_;

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
 
   /* 
   inline double LocalLamellarOrderingExternal::energy(double dPerp, double dParallel, int type) const
   {
      double perpLength, parallelLength, q, clipParameter, arg, clipcos, parallelFactor;
      Vector lengths;
      lengths = boundaryPtr_->lengths();
      perpLength = lengths[perpDirection_];
      q = (2.0*M_PI*periodicity_)/perpLength;
      clipParameter   = 1.0/(q*width_*perpLength);
      arg = q*dPerp;
      clipcos = clipParameter*cos(arg);
      parallelLength = lengths[parallelDirection_];
      parallelFactor = tanh( 2.0*M_PI*((dParallel - (((1.0-fraction_)/2.0)*parallelLength))/(fraction_*parallelLength)) );
      parallelFactor *= tanh( 2.0*M_PI*(((((1.0+fraction_)/2.0)*parallelLength) - dParallel)/(fraction_*parallelLength)) );
      parallelFactor += 1.0;
      parallelFactor /= 2.0;
      if (parallelFactor < 0.0)
         UTIL_THROW("Negative base factor");
      return prefactor_[type]*externalParameter_*tanh(clipcos)*parallelFactor;
   }
   */

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline 
   double LocalLamellarOrderingExternal::energy(const Vector& position, int type) const
   {
      double dPerp, dParallel, perpLength, parallelLength, q, clipParameter, arg, clipcos;
      double parallelFactor1, parallelFactor2, parallelFactor;
      dPerp = position[perpDirection_];
      dParallel = position[parallelDirection_];

      Vector lengths;
      lengths = boundaryPtr_->lengths();
      perpLength = lengths[perpDirection_];
      parallelLength = lengths[parallelDirection_];

      q = (2.0*M_PI*periodicity_)/perpLength;
      clipParameter   = 1.0/(q*width_*perpLength);
      arg = q*dPerp;
      clipcos = clipParameter*cos(arg);
      parallelFactor1 = tanh( 2.0*M_PI*((dParallel - (((1.0-fraction_)/2.0)*parallelLength))/(fraction_*parallelLength)) );
      parallelFactor2 = tanh( 2.0*M_PI*(((((1+fraction_)/2.0)*parallelLength) - dParallel)/(fraction_*parallelLength)) );
      parallelFactor = 1.0 + (parallelFactor1*parallelFactor2);
      parallelFactor /= 2.0;

      if (parallelFactor < 0.0)
         UTIL_THROW("Negative base factor");

      return prefactor_[type]*externalParameter_*tanh(clipcos)*parallelFactor;
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void LocalLamellarOrderingExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      double dPerp, dParallel, perpLength, parallelLength, q, clipParameter, arg, clipcos, tanH, sechSq;
      double parallelFactor1, parallelFactor2, parallelFactor, forceParallelSech1, forceParallelSech2;
      dPerp = position[perpDirection_];
      dParallel = position[parallelDirection_];
      force.zero();

      Vector lengths;
      lengths = boundaryPtr_->lengths();
      perpLength = lengths[perpDirection_];
      parallelLength = lengths[parallelDirection_];

      q = (2.0*M_PI*periodicity_)/perpLength;
      clipParameter  = 1.0/(q*width_*perpLength);
      arg = q*dPerp;
      clipcos = clipParameter*cos(arg);
      tanH = tanh(clipcos);
      sechSq = (1.0 - tanH*tanH);
      parallelFactor1 = tanh( 2.0*M_PI*((dParallel - (((1.0-fraction_)/2.0)*parallelLength))/(fraction_*parallelLength)) );
      parallelFactor2 = tanh( 2.0*M_PI*(((((1+fraction_)/2.0)*parallelLength) - dParallel)/(fraction_*parallelLength)) );
      parallelFactor = 1.0 + (parallelFactor1*parallelFactor2);
      parallelFactor /= 2.0;

      if (parallelFactor < 0.0) {
         UTIL_THROW("Negative base factor");
      } else {
         forceParallelSech1 = (1.0 - (parallelFactor1*parallelFactor1));
         forceParallelSech1 *= 1.0*M_PI*(1/(fraction_*parallelLength))*parallelFactor2;
         forceParallelSech2 = (1.0 - (parallelFactor2*parallelFactor2));
         forceParallelSech2 *= -1.0*M_PI*(1/(fraction_*parallelLength))*parallelFactor1;
         force[parallelDirection_] = -1.0*prefactor_[type]*externalParameter_*tanh(clipcos);
         force[parallelDirection_] *= (forceParallelSech1 + forceParallelSech2);
         force[perpDirection_] = prefactor_[type]*externalParameter_*sechSq*clipParameter*sin(arg)*q*parallelFactor;
      }

   }
 
}
#endif
