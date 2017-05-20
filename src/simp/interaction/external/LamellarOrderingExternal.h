#ifndef SIMP_LAMELLAR_ORDERING_EXTERNAL_H
#define SIMP_LAMELLAR_ORDERING_EXTERNAL_H

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
   * \ingroup Simp_Interaction_External_Module
   */
   class LamellarOrderingExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      LamellarOrderingExternal();

      /**
      * Copy constructor.
      */
      LamellarOrderingExternal(const LamellarOrderingExternal& other);

      /**
      * Assignment.
      */
      LamellarOrderingExternal& operator = (const LamellarOrderingExternal& other);

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
      * Return name string "LamellarOrderingExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 2;

      /// Index representing the direction perpendicular to the lamellae.
      int perpDirection_;

      /// Interfcial width in lamellar phase.
      double interfaceWidth_;

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
 
    
   inline double LamellarOrderingExternal::energy(double d, int type) const
   {
      double perpLength, q, clipParameter, arg, clipcos;
      Vector lengths;
      lengths = boundaryPtr_->lengths();
      perpLength = lengths[perpDirection_];

      q = (2.0*M_PI*periodicity_)/perpLength;
      clipParameter   = 1.0/(q*interfaceWidth_*perpLength);
      arg = q*d;
      clipcos = clipParameter*cos(arg);
      
      return prefactor_[type]*externalParameter_*tanh(clipcos);
   }

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline 
   double LamellarOrderingExternal::energy(const Vector& position, int type) const
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
   inline double LamellarOrderingExternal::forceScalar(double d, int type) const
   {
      double perpLength, q, clipParameter, arg, clipcos, tanH, sechSq;
      Vector lengths;
      lengths = boundaryPtr_->lengths();
      perpLength = lengths[perpDirection_];

      q = (2.0*M_PI*periodicity_)/perpLength;
      clipParameter = 1.0/(q*interfaceWidth_*perpLength);
      arg = q*d;
      clipcos = clipParameter*cos(arg);
      tanH = tanh(clipcos);
      sechSq = (1.0 - tanH*tanH);
      return prefactor_[type]*externalParameter_*sechSq*clipParameter*sin(arg)*q;
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void LamellarOrderingExternal::getForce(const Vector& position, int type, 
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
