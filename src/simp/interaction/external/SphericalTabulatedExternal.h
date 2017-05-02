#ifndef SIMP_SPHERICAL_TABULATED_EXTERNAL_H
#define SIMP_SPHERICAL_TABULATED_EXTERNAL_H

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
#include <string>

namespace Simp
{

   using namespace Util;

   /**
   * The potential is independent on the type of atoms.
   * 
   * \ingroup Simp_Interaction_External_Module
   */
   class SphericalTabulatedExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      SphericalTabulatedExternal();

      /**
      * Copy constructor.
      */
      SphericalTabulatedExternal(const SphericalTabulatedExternal& other);

      /**
      * Assignment.
      */
      SphericalTabulatedExternal& 
      operator = (const SphericalTabulatedExternal& other);

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
      * Return name string for this interaction class.
      */
      std::string className() const;

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
 
      //@}

   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 6;

      /// Filename for tabulated potentials
      std::string filename_;

      /// Cutoff distance.
      double rMax_;

      /// Square of cutoff distance
      double rMaxSq_;

      /// Grid step size
      double dr_;
  
      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of grid points (including endpoints).
      int nr_; 

      /// Number of possible atom types.
      int nAtomType_; 

      /// Tabulated potentials (one per atom type).
      DArray< DArray<double> > potentials_;

      /// Are all parameters and pointers initialized?
      bool isInitialized_;

   };
  
   // inline methods 
 
   /* 
   * Calculate external potential energy for a single atom.
   */
   inline 
   double SphericalTabulatedExternal::energy(const Vector& position, int type) const
   {
      Vector origin = boundaryPtr_->lengths();
      origin *= 0.5;
      double rSq = boundaryPtr_->distanceSq(position, origin);
      if (rSq < rMaxSq_) {
         double r = sqrt(rSq);
         int k = (int) (r/dr_);
         double fp = r - k*dr_;
         UTIL_CHECK(k >= 0);
         UTIL_CHECK(k < nr_ -1);
         UTIL_CHECK(fp >= 0.0);
         UTIL_CHECK(fp <= 1.0);
         double fm = 1.0 - fp;
         return fm*potentials_[type][k] + fp*potentials_[type][k+1]; 
      } else {
         return 0.0;
      }
   }
   
   /* 
   * Calculate external force for a single atom.
   */
   inline void 
   SphericalTabulatedExternal::getForce(const Vector& position, int type, 
                                        Vector& force) 
   const
   {
      Vector origin = boundaryPtr_->lengths();
      origin *= 0.5;
      double rSq = boundaryPtr_->distanceSq(position, origin, force);
      if (rSq < rMaxSq_) {
         double r = sqrt(rSq);
         int k = (int) (r/dr_);
         UTIL_CHECK(k >= 0);
         UTIL_CHECK(k < nr_ -1);
         double fs = potentials_[type][k] - potentials_[type][k+1];
         fs /= dr_;
         force *= fs/r;
      } else {
         for (int i = 0; i < Dimension; ++i) {
            force[i] = 0.0;
         }
      }
   }
 
}
#endif
