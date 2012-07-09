#ifndef DDMD_PAIR_POTENTIAL_H
#define DDMD_PAIR_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
#include <ddMd/neighbor/CellList.h>     // member
#include <ddMd/neighbor/PairList.h>     // member
#include <ddMd/util/DdTimer.h>          // member
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
   */
   class PairPotential : public ParamComposite
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
      * Required if this was created with default constructor PairPotential().
      * The associate method must be called before setParam. This method sets 
      * the skin and cutoff length parameters, and allocates the internal
      * CellList and a PairList. The CellList allocates memory for the number
      * of cells required for the specified upper and lower coordinate bounds, 
      * with the specified cutoff. 
      *
      * \param lower        lower dimension of the largest expected domain
      * \param upper        upper dimension of the largest expected domain
      * \param skin         pair list skin length
      * \param pairCapacity maximum number of pairs per processor
      */
      void initialize(const Vector& lower, const Vector& upper,
                      double skin, int pairCapacity);

      /**
      * Set flag to identify if reverse communication is enabled.
      *
      * \param reverseUpdateFlag true if reverse communication is enabled.
      */
      void setForceCommFlag(bool reverseUpdateFlag);

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
      * Return pair energy for a single pair.
      * 
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual 
      double energy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
      */
      virtual 
      double forceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add pair forces to atom forces.
      */
      virtual void addForces() = 0;

      /**
      * Calculate total pair potential on all processors.
      *
      * This method must be called on all processors. The result is
      * stored on the master processor, and may be retrieved by 
      * calling energy() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy() = 0;
      #endif

      /**
      * Return the total pair potential, on all processors.
      *
      * This method should only be called on the master (rank 0) processor,
      * after a previous call to computeEnergy.
      */
      virtual double energy() = 0;

      #if 0
      /**
      * Compute total bond pressure.
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z bond pressure components.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute bond stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;
      #endif

      //@}
      /// \name Pair and Cell Lists.
      //@{

      /**
      * Build the cell and pair lists. 
      *
      * Use with objects created with PairPotential(Simulation&). Makes a
      * CellList grid for the region defined by the associated Domain. 
      */
      void findNeighbors();

      /**
      * Build the cell and pair lists. 
      *
      * \param lower Vector of lower coordinate bounds for this processor.
      * \param upper Vector of upper coordinate bounds for this processor.
      */
      void findNeighbors(const Vector& lower, const Vector& upper);

      /**
      * Compute twice the number of pairs within force cutoff, on all processors.
      *  
      * This method must be called on all processors.
      *
      * This method add +1 for each atom in such a pair, or +2 for a complete
      * pair. The method for distributing pairs among processors is the same
      * as is used to distribute the energy, and depends on the value of 
      * reverseUpdateFlag().
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
      * Is reverse communication enabled?
      */
      bool reverseUpdateFlag() const;

      /**
      * Return integer id for algorithm (0=PAIR, 1=CELL, 2=NSQ)
      */
      int methodId() const;

      /**
      * Return internal timer by reference
      */
      DdTimer& timer();

      /**
      * Enumeration of time stamp identifiers.
      */
      enum TimeId {START, BUILD_CELL_LIST, BUILD_PAIR_LIST, 
                   FORCES, NTime};

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
      * Read parameters and allocate memory for PairList.
      *
      * Use iff this object was instantiated with PairPotential(Simulation&).
      *
      * \param in input parameter stream.
      */
      void readPairListParam(std::istream& in);

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
      * Stamp timer with int/enum TimeId.
      */
      void stamp(unsigned int timeId);

   private:

      // Pointer to associated Domain object.
      Domain* domainPtr_;

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* storagePtr_;

      /// Timer
      DdTimer timer_;

      /// Index for method used to calculate forces / energies.
      int methodId_;

      /// Number of pairs within specified cutoff.
      int nPair_;

      /// Is reverse communication enabled?
      bool reverseUpdateFlag_;

      // Private methods used to compute number of pairs
      int nPairList(double cutoffSq);
      int nPairCell(double cutoffSq);
      int nPairNSq(double cutoffSq);

   };

   inline CellList& PairPotential::cellList()
   {  return cellList_; }

   inline PairList& PairPotential::pairList()
   {  return pairList_; }

   inline double PairPotential::skin() const
   {  return skin_; }

   inline double PairPotential::cutoff() const
   {  return cutoff_; }

   inline bool PairPotential::reverseUpdateFlag() const
   {  return reverseUpdateFlag_; }

   inline Boundary& PairPotential::boundary() 
   {  return *boundaryPtr_; }

   inline Domain& PairPotential::domain()
   {  return *domainPtr_; }

   inline AtomStorage& PairPotential::storage()
   {  return *storagePtr_; }

   inline DdTimer& PairPotential::timer()
   {  return timer_; }

   inline void PairPotential::stamp(unsigned int timeId) 
   {
      #ifdef DDMD_PAIR_POTENTIAL_TIMER
      timer_.stamp(timeId);
      #endif
   }

   inline void PairPotential::setMethodId(int methodId)
   {  methodId_ = methodId; }

   inline int PairPotential::methodId() const
   {  return methodId_; }

}
#endif
