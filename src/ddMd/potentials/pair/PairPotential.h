#ifndef DDMD_PAIR_POTENTIAL_H
#define DDMD_PAIR_POTENTIAL_H

#include <ddMd/potentials/Potential.h>  // base class
#include <ddMd/neighbor/CellList.h>     // member
#include <ddMd/neighbor/PairList.h>     // member
#include <ddMd/misc/DdTimer.h>          // member
#include <util/boundary/Boundary.h>     // member (typedef)
#include <util/global.h>

#include <iostream>

#define DDMD_PAIR_POTENTIAL_TIMER

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{  class Tensor; }

namespace DdMd
{

   class Simulation;
   class AtomStorage;
   class Domain;
   using namespace Util;

   /**
   * A PairPotential calculates forces for a parent Simulation.
   *
   * A PairPotential has a private CellList and PairList which it
   * uses to calculate nonbonded pair forces. 
   *
   * \ingroup DdMd_Pair_Module
   */
   class PairPotential : public Potential
   {

   public:

      /**
      * Default constructor (for unit testing).
      */
      PairPotential();

      /**
      * Constructor.
      */
      PairPotential(Simulation& simulation);

      /**
      * Associate with related objects.
      *
      * Call iff object instantiated with default constructor.
      *
      * \param domain   associated Domain object.
      * \param boundary associated Boundary object.
      * \param storage  associated AtomStorage object.
      */
      void associate(Domain& domain, Boundary& boundary, AtomStorage& storage);

      /**
      * Destructor.
      */
      virtual ~PairPotential();

      /**
      * Set parameters and allocate memory.
      *
      * This method sets values for the same members as readParameters.
      * Either method must be called after associate. This method sets the
      * skin and cutoff length parameters, and allocates memory for the 
      * internal CellList and a PairList. It uses the maximum boundary and
      * the cutoff to calculate the how many cells to allocate.
      *
      * \param maxBoundary  largest expected Boundary (used for allocation)
      * \param skin         pair list skin length
      * \param pairCapacity maximum number of pairs per processor
      */
      void 
      initialize(const Boundary& maxBoundary, double skin, int pairCapacity);

      /**
      * Read parameters and allocate memory for PairList.
      *
      * Use iff this was instantiated with PairPotential(Simulation&).
      *
      * \param in input parameter stream.
      */
      void readParameters(std::istream& in);

      /**
      * Load parameters for PairList from archive, and allocate memory.
      *
      * Use iff this was instantiated with PairPotential(Simulation&).
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save parameters for PairList to output archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Set id to specify algorithm for energy, force calculations.
      *
      * \param methodId algorithm id: 0=pair list, 1=cell list, 2=N^2 loop.
      */
      void setMethodId(int methodId);

      /// \name Interaction interface
      //@{

      /**
      * Set the maximum number of atom types.
      */
      virtual void setNAtomType(int nAtomType) = 0;
 
      /**
      * Return energy for a single pair.
      * 
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual 
      double pairEnergy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
      */
      virtual 
      double pairForceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  new value of parameter
      */
      virtual void set(std::string name, int i, int j, double value) = 0;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      virtual double get(std::string name, int i, int j) const = 0;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name Pair and Cell Lists.
      //@{

      /**
      * Build a cell list.
      *
      * This method adds all local and ghost atoms to a cell list that is
      * used in buildPairList() to build a Verlet pair list. On entry to
      * this method, atomic positions must Cartesian if UTIL_ORTHOGONAL
      * is true and generalized coordinates if UTIL_ORTHOGONAL is false.
      * If UTIL_ORTHOGONAL is false, this method should be followed by a 
      * call to AtomStorage::transformGenToCart(Boundary& ) to transform
      * all positions back to Cartesian coordinates.
      */
      void buildCellList();

      /**
      * Build the Verlet Pair list.
      *
      * Preconditions:
      *  - Cell List must have been built, by PairPotential::buildCellList.
      *  - Atomic positions must be expressed in Cartesian coordinates.
      */
      void buildPairList();

      /**
      * Compute pair energies on all processors.
      *
      * This method must be called on all processors. The result is
      * stored on the master processor, and may be retrieved by 
      * calling pairEnergies() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computePairEnergies(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computePairEnergies() = 0;
      #endif

      /**
      * Return total pair energies, from all processors.
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computePairEnergies.
      */
      Util::DMatrix<double> pairEnergies() const;

      /**
      * Mark pair energy as unknown (nullify).
      */
      void unsetPairEnergies();

      /**
      * Compute twice the number of pairs within the force cutoff.
      *  
      * This method must be called on all processors.
      *
      * The method for distributing pairs among processors is the same
      * as is used to distribute the energy, and depends on the value
      * of reverseUpdateFlag().
      */
      #ifdef UTIL_MPI
      void computeNPair(MPI::Intracomm& communicator);
      #else
      void computeNPair();
      #endif

      /**
      * Return twice the number of pairs within the specified force cutoff.
      * 
      * This method should only be called on the rank 0 processor. The
      * return value is computed by a previous call to computeNPair.
      */
      int nPair() const;

      /**
      * Get the CellList by reference.
      */
      CellList& cellList();

      /**
      * Get the PairList by reference.
      */
      PairList& pairList();

      //@}

      /**
      * Get value of the pair list skin.
      */
      double skin() const;

      /**
      * Get value of the pair list cutoff (maxPairCutoff + skin).
      */
      double cutoff() const;

      /**
      * Return integer id for algorithm (0=PAIR, 1=CELL, 2=NSQ)
      */
      int methodId() const;

   protected:

      // CellList to construct PairList or calculate nonbonded pair forces.
      CellList cellList_;

      // Verlet pair list, to calculate nonbonded pair forces.
      PairList pairList_;

      // Boundary used to allocate space for the cell list.
      Boundary maxBoundary_;

      // Difference between pairlist cutoff and pair potential cutoff. 
      double skin_;

      // Minimum cell size = pair potential cutoff + skin.
      double cutoff_;

      // Maximum number of nonbonded pairs in pair list. 
      int pairCapacity_;

      /**
      * Get the PairList by const reference.
      */
      Boundary& boundary();

      /**
      * Get the PairList by const reference.
      */
      Domain& domain();

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& storage();

      /**
      * Set values for pair energies.
      */
      void setPairEnergies(DMatrix<double> pairEnergies);

   private:

      // Pointer to associated Domain object.
      Domain* domainPtr_;

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

      /// Index for method used to calculate forces / energies.
      int methodId_;

      /// Number of pairs within specified cutoff.
      int nPair_;

      /// Pair energies.
      Setable< DMatrix<double> > pairEnergies_;

      // Private methods used to compute number of pairs
      int nPairList(double cutoffSq);
      int nPairCell(double cutoffSq);
      int nPairNSq(double cutoffSq);

      /*
      * Allocate memory for the cell list and pair list.
      */
      void allocate();

   };

   inline CellList& PairPotential::cellList()
   {  return cellList_; }

   inline PairList& PairPotential::pairList()
   {  return pairList_; }

   inline double PairPotential::skin() const
   {  return skin_; }

   inline double PairPotential::cutoff() const
   {  return cutoff_; }

   inline Boundary& PairPotential::boundary() 
   {  return *boundaryPtr_; }

   inline Domain& PairPotential::domain()
   {  return *domainPtr_; }

   inline AtomStorage& PairPotential::storage()
   {  return *storagePtr_; }

   inline void PairPotential::setMethodId(int methodId)
   {  methodId_ = methodId; }

   inline int PairPotential::methodId() const
   {  return methodId_; }

}
#endif
