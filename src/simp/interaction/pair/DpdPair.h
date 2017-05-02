#ifndef SIMP_DPD_PAIR_H
#define SIMP_DPD_PAIR_H

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
   * Soft pair potential used in dissipative particle dynamics (DPD) simulations 
   * of Groot, Warren et al. 
   *
   * \sa \ref simp_interaction_pair_DpdPair_page "Parameter file format"
   * \sa \ref simp_interaction_pair_interface_page
   * \sa \ref simp_interaction_pair_page
   * 
   * \ingroup Simp_Interaction_Pair_Module
   */
   class DpdPair : public ParamComposite 
   {
   
   public:
   
      /**
      * Constructor.
      */
      DpdPair();

      /**
      * Copy constructor.
      */
      DpdPair(const DpdPair& other);

      /**
      * Assignment.
      */
      DpdPair& operator = (const DpdPair& other);

      /// \name Mutators
      //@{ 

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Read epsilon and sigma, initialize other variables.
      *
      * \pre nAtomType must be set, by calling setNAtomType().
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
      /// \name Accessors (required)
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
      * Precondition: The square separation rsq must be less than cutoffSq.
      * If rsq > cutoffSq, the return value is undefined (i.e., wrong).
      * Usage: Test for rsq < cutoffSq before calling this function
      * \code
      * if (rsq < interaction.cutoffSq(i, j)) {
      *    f = forceOverR(rsq, i, j);
      *    .....
      * }
      * \endcode
      *
      * \param rsq square of distance between particles
      * \param i type of particle 1
      * \param j type of particle 2
      * \return  force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;
   
      /**
      * Get square of cutoff distance for specific type pair.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    cutoffSq_[i][j]
      */
      double cutoffSq(int i, int j) const;
 
      /**
      * Get maximum of pair cutoff distance, for all atom type pairs.
      */
      double maxPairCutoff() const;
 
      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      */
      double get(std::string name, int i, int j) const;

      //@}
      /// \name Accessors (extra)
      //@{

      /**
      * Get interaction energy for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double epsilon(int i, int j) const;
 
      /**
      * Get range for a specific pair of Atom types.
      *
      * \param i   atom type index 1
      * \param j   atom type index 2
      * \return    sigma_[i][j]
      */
      double sigma(int i, int j) const;
 
      //@}

   private:
   
      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxAtomType = 4;
   
      // Parameters for different types of particle pairs
      double epsilon_[MaxAtomType][MaxAtomType];   ///< energy parameter
      double sigma_[MaxAtomType][MaxAtomType];     ///< range parameter
      double sigmaSq_[MaxAtomType][MaxAtomType];   ///< square of sigma[][].
      double ce_[MaxAtomType][MaxAtomType];        ///< energy prefactor
      double cf_[MaxAtomType][MaxAtomType];        ///< force prefactor
 
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
   inline double DpdPair::energy(double rsq, int i, int j) const 
   {
      double dr;
      if (rsq < sigmaSq_[i][j]) {
         dr = sqrt(rsq) - sigma_[i][j];
         return ce_[i][j]*dr*dr;
      } else {
         return 0.0;
      }
   }
   
   /* 
   * Calculate force/distance for a pair as function of squared distance.
   */
   inline double DpdPair::forceOverR(double rsq, int i, int j) const
   {
      if (rsq < sigmaSq_[i][j]) {
         return cf_[i][j]*(sigma_[i][j]/sqrt(rsq) - 1.0);
      } else {
         return 0.0;
      }
   }

   /* 
   * Calculate force/distance for a pair as function of squared distance.
   */
   inline double DpdPair::cutoffSq(int i, int j) const
   {  return sigmaSq_[i][j]; }

}
#endif
