#ifndef INTER_LJ_PAIR_H
#define INTER_LJ_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

namespace Inter
{

   using namespace Util;

   /**
   * A cutoff Lennard-Jones nonbonded pair potential.
   *
   * Member variables include interaction parameters. Member functions evaluate
   * energy and force for an individual pair of nonbonded interacting particles.
   * 
   * \ingroup Inter_Pair_Module
   */
   class LJPair : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      LJPair();

      /**
      * Copy constructor.
      */
      LJPair(const LJPair& other);

      /**
      * Assignment.
      */
      LJPair& operator = (const LJPair& other);

      /// \name Mutators
      //@{ 

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Read epsilon, sigma, and cutoff, and initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
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

      /**
      * Set LJ interaction energy for a specific pair of Atom types.
      *
      * \param i        type of Atom 1
      * \param j        type of Atom 2
      * \param epsilon  LJ energy parameter
      */
      void setEpsilon(int i, int j, double epsilon);
 
      /**
      * Get LJ range for a specific pair of Atom types.
      *
      * \param i      atom type index 1
      * \param j      atom type index 2
      * \param sigma  LJ range parameter
      */
      void setSigma(int i, int j, double sigma);

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value);

      //@}
      /// \name Accessors
      //@{ 

      /**
      * Returns interaction energy for a single pair of particles. 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    pair interaction energy
      */
      double energy(double rsq, int i, int j) const;
  
      /**
      * Returns ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of particles.
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;
  
      /**
      * Get maximum of pair cutoff distance, for all atom type pairs.
      */
      double maxPairCutoff() const;

      /**
      * Get LJ interaction energy for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double epsilon(int i, int j) const;
 
      /**
      * Get LJ range for a specific pair of Atom types.
      *
      * \param i   atom type index 1
      * \param j   atom type index 2
      * \return    sigma_[i][j]
      */
      double sigma(int i, int j) const;
 
      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      */
      double get(std::string name, int i, int j) const;

      //@}

   private:
   
      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxAtomType = 2;
   
      // Lennard-Jones parameters for different types of particle pairs
      double epsilon_[MaxAtomType][MaxAtomType];   ///< LJ interaction energies.
      double sigma_[MaxAtomType][MaxAtomType];     ///< LJ range parameters.
      double sigmaSq_[MaxAtomType][MaxAtomType];   ///< square of sigma[][].
      double cutoff_[MaxAtomType][MaxAtomType];    ///< LJ cutoff distance.
      double cutoffSq_[MaxAtomType][MaxAtomType];  ///< square of cutoff[][].
      double ljShift_[MaxAtomType][MaxAtomType];   ///< shift in LJ potential.
 
      /**
      * Maximum pair potential cutoff radius, for all monomer type pairs.
      *
      * Used in construction of a cell list or Verlet pair list.
      */
      double maxPairCutoff_;

      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 
 
   /* 
   * Calculate interaction energy for a pair, as function of squared distance.
   */
   inline double LJPair::energy(double rsq, int i, int j) const 
   {
      double r6i;
      if (rsq < cutoffSq_[i][j]) {
         if (rsq < 0.6*sigmaSq_[i][j]) {
             rsq = 0.6*sigmaSq_[i][j];
         }
         r6i = sigmaSq_[i][j]/rsq;
         r6i = r6i*r6i*r6i;
         return 4.0*epsilon_[i][j]*(r6i*r6i - r6i) + ljShift_[i][j];
      } else {
         return 0.0;
      }
   }
  
   /* 
   * Calculate force/distance for a pair as function of squared distance.
   */
   inline double LJPair::forceOverR(double rsq, int i, int j) const
   {
      double r2i, r6i;
      if ( rsq < cutoffSq_[i][j] ) {
         if ( rsq < 0.6*sigmaSq_[i][j] ) {
             return 0.0;
         }
         r2i = 1.0/rsq;
         r6i = sigmaSq_[i][j]*r2i;
         r6i = r6i*r6i*r6i;
         return 24.0*epsilon_[i][j]*(2.0*r6i*r6i - r6i)*r2i;
      } else {
         return 0.0;
      }
   }

}
#endif
