#ifndef SIMP_LJ_PAIR_H
#define SIMP_LJ_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

namespace Simp
{

   using namespace Util;

   /**
   * A cutoff, shifted Lennard-Jones nonbonded pair interaction.
   *
   * This class defines a Lennard-Jones potential that is cutoff and
   * shifted so that the potential is zero at the cutoff distance.
   *
   * \sa \ref simp_interaction_pair_LJPair_page "Parameter file format"
   * \sa \ref simp_interaction_pair_interface_page
   * \sa \ref simp_interaction_pair_page
   * 
   * \ingroup Simp_Interaction_Pair_Module
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
      *
      * \param other LJPair to be copied
      */
      LJPair(const LJPair& other);

      /**
      * Assignment.
      *
      * \param other LJPair to be assigned.
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
      * \param in  input parameter stream 
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
      * \param i      atom type index (1st atom)
      * \param j      atom type index (2nd atom)
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
      * Returns interaction energy for a single pair of atoms. 
      *
      * \param rsq square of distance between atoms
      * \param i   type of atom 1
      * \param j   type of atom 2
      * \return    pair interaction energy
      */
      double energy(double rsq, int i, int j) const;
  
      /**
      * Returns ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of atoms.
      *
      * Precondition: The distance squared rsq must be less than cutoffSq.
      * If rsq > cutoffSq, the return value is undefined (i.e., invalid).
      * Usage: Test for rsq < cutoffSq before calling this function
      * \code
      * if (rsq < interaction.cutoffSq(i, j)) {
      *    f = forceOverR(rsq, i, j);
      *    .....
      * }
      * \endcode
      *
      * \param rsq square of distance between atoms
      * \param i   type of atom 1
      * \param j   type of atom 2
      * \return    force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;
  
      /**
      * Get square of cutoff distance for specific type pair.
      *
      * \param i   type of atom 1
      * \param j   type of atom 2
      * \return    cutoffSq_[i][j]
      */
      double cutoffSq(int i, int j) const;
 
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

   protected:
   
      /// Maximum allowed value for nAtomType (# of atom types)
      static const int MaxAtomType = 4;
   
      // Lennard-Jones parameters for different types of atom pairs
      double epsilon_[MaxAtomType][MaxAtomType];  ///< LJ interaction energies.
      double sigma_[MaxAtomType][MaxAtomType];    ///< LJ range parameters.
      double sigmaSq_[MaxAtomType][MaxAtomType];  ///< square of sigma[][].
      double cutoff_[MaxAtomType][MaxAtomType];   ///< LJ cutoff distance.
      double cutoffSq_[MaxAtomType][MaxAtomType]; ///< square of cutoff[][].
      double ljShift_[MaxAtomType][MaxAtomType];  ///< shift in LJ potential.
      double eps48_[MaxAtomType][MaxAtomType];    ///< 48*epsilon
 
      /**
      * Maximum pair potential cutoff radius, for all monomer type pairs.
      *
      * Used in construction of a cell list or Verlet pair list.
      */
      double maxPairCutoff_;

      /**
      * Total number of atom types.
      */
      int    nAtomType_; 

      /**
      * Was this object initialized by calling (read|load)Parameters ?
      */
      bool  isInitialized_;

   };
  
   // Inline methods 
 
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
      r2i = 1.0/rsq;
      r6i = sigmaSq_[i][j]*r2i;
      r6i = r6i*r6i*r6i;
      return eps48_[i][j]*(r6i - 0.5)*r6i*r2i;
   }

   /* 
   * Return cutoff parameter for a specific atom type pair.
   */
   inline double LJPair::cutoffSq(int i, int j) const
   {  return cutoffSq_[i][j]; }

}
#endif
