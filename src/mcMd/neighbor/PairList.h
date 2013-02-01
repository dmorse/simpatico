#ifndef MCMD_PAIR_LIST_H
#define MCMD_PAIR_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CellList.h"
#include <util/boundary/Boundary.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>
#include <util/space/Vector.h>

class PairListTest;

namespace McMd
{

   using namespace Util;

   class Atom;
   class PairIterator;
   
   /**
   * A Verlet neighbor list.
   *
   * A PairList (or Verlet list) is a list of neighboring pairs of Atoms that 
   * are separated by a distance less than a specified cutoff. The cutoff for
   * the Verlet list is the sum of a potential cutoff, which is passed as a 
   * parameter to allocate(), and a "skin", which is read by readParameters().
   *
   * After a PairList is constructed, the allocate() method must be called to
   * allocate memory for both the data structures required to store the PairList
   * and for a private CellList object that is used to construct the PairList.
   *
   * To build or rebuild a PairList, after it has been allocated, one must 
   * invoke clear() to clear the private CellList, call addAtom(Atom&) once for
   * every atom in the System, within a loop over atoms,  and then call build() 
   * to actually construct the PairList. The addAtom() method adds an Atom to
   * an internal CellList, and then build() uses the CellList to build a list
   * of pairs. Atomic positions must not be changed during this process.
   *
   * A PairIterator object must be used to iterate over all of the pairs in a 
   * completed PairList (see documentation for PairIterator for usage).
   *
   * A completed PairList is guaranteed to remain valid as long as no Atom 
   * in the system has moved more than a distance skin/2 from its position 
   * when the PairList was last built. The method isCurrent() determines if 
   * this criterion is satisfied, by comparing current atomic positions to
   * positions that were stored when the PairList was last built. In an MD
   * simulation, isCurrent() should be called after every time step, and 
   * the PairList should be rebuilt if its return value is false.
   * 
   *
   * \ingroup McMd_Neighbor_Module 
   */
   class PairList : public ParamComposite 
   {

      // Null value for any non-negative index
      static const int NullIndex = -1;   

   public:
  
      /**
      * Default constructor.
      */
      PairList();
   
      /**
      * Destructor. 
      *
      * Deallocates arrays atom1Ptrs_, atom2Ptrs_, first_, and oldPositions_.
      */
      virtual ~PairList();
   
      /**
      * Read atomCapacity and pairCapacity (maximum numbers of atoms and pairs).
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

      /// \name Mutators
      //@{

      /**
      * Allocate memory and initialize.
      *
      * Allocate all memory required by the PairList and the CellList,
      * and initialize the grid of Cell objects. 
      *
      * Precondition: readParameters() must be invoked before allocate(),
      * so that values of atomCapacity, pairCapacity, skin are known.
      * The allocate() method can only be called once.
      *
      * \param atomIdEnd       maximum allowed atom Id, plus 1
      * \param boundary        Boundary object with maximum dimensions
      * \param potentialCutoff Range of pair potential
      */
      void 
      allocate(int atomIdEnd, const Boundary &boundary, double potentialCutoff);
  
      /**
      * Make the grid of cells for the internal Cell List.
      *
      * Precondition: This PairList must be allocated.
      */ 
      void makeGrid(const Boundary &boundary);

      /**
      * Clear the PairList and CellList.
      *
      * Precondition: This PairList must be allocated.
      */ 
      void clear();

      /**
      * Add an Atom to the CellList.
      *
      * After the clear() method is invoked, every Atom in a System must 
      * be passed to this method sequentially, within a loop over Atoms.
      * The Atom positions must not change within this loop. 
      *
      * \param atom Atom to be added.
      */ 
      void addAtom(Atom& atom);

      /**
      * Use a complete CellList to build a new PairList.
      *
      * Preconditions: 
      *   - This PairList must be allocated
      *   - All atoms must have been added to the CellList.
      *   - Atomic positions must not have changed since they were added.
      *
      * \param boundary current Boundary object. 
      */
      void build(const Boundary& boundary);

      /**
      * Initialize a PairIterator.
      *
      * \param iterator a PairList, initialized on output
      */
      void begin(PairIterator &iterator) const;
 
      //@}
      /// \name Accessors (miscellaneous)
      //@{ 
 
      /**
      * Get the number of pairs currently in the PairList.
      */
      int nAtom() const;

      /**
      * Get the number of pairs currently in the PairList.
      */
      int nPair() const;

      /**
      * Has memory been allocated for this PairList?
      */
      bool isAllocated() const;
   
      /**
      * Returns true if PairList is current, false otherwise.
      *
      * The PairList is considered current if no atom has moved a distance
      * skin/2 since the PairList was built, and obsolete otherwise. The
      * parameter skin is the difference between the Verlet radius that 
      * was used to construct the PairList and the maximum range of any 
      * nonbonded potential (i.e., the maxCutoff member of a PairPotential 
      * object).
      *
      * \param  boundary  Boundary object containing simulation cell dimensions
      * \return true if the Pairlist is valid, false if it is outdated
      */
      bool isCurrent(const Boundary &boundary) const;

      /// Return true if valid, or throw Exception.
      bool isValid() const;
  
      //@}
      /// \name Statistics
      //@{
 
      /**
      * Get the maximum number of atoms encountered thus far.
      */
      int maxNAtom() const;

      /**
      * Get the maximum number of pairs encountered thus far.
      */
      int maxNPair() const;

      /**
      * Return number of times the PairList has been built thus far.
      */
      int buildCounter() const;

      /**
      * Clear statistical accumulators (maxNAtom, maxNPair, buildCounter).
      */
      void clearStatistics();

      //@}

   private:
  
      /// Private CellList, used to create PairList.
      CellList cellList_;

      /// Array of pointers to 1st (or primary) atom in each pair.
      DArray<Atom*>  atom1Ptrs_;  

      /// Array of pointers to neighbor (or secondary) atom in each pair.
      DArray<Atom*>  atom2Ptrs_;  

      /// Array of indices in atom2Ptrs_ of first neighbor of an Atom.
      DArray<int>    first_; 

      /// Array of old atom positions.
      DArray<Vector> oldPositions_;

      /// Extra distance to add to pair potential cutoff.
      double skin_;
   
      /// Pair list cutoff radius (pair potential cutoff + skin_).
      double cutoff_;
   
      /// Maximum number of atoms (dimension of atom1Ptrs_).
      int    atomCapacity_;     
   
      /// Maximum number of distinct pairs (dimension of atom2Ptrs_).
      int    pairCapacity_;     
   
      /// Number of primary atoms in atom1Ptrs_.
      int    nAtom1_;      

      /// Number of secondary atoms in atom2Ptrs_, or number of pairs.
      int    nAtom2_; 
   
      /// Total number of atoms in the PairList.
      int    nAtom_;     
   
      /// Index one less than the first element of atom1Ptrs_ for atoms 
      /// with no neighbors.
      int    tList1_;     

      /// Maximum value of nAtom_ since instantiation.
      int    maxNAtom_;     
   
      /// Maximum value of nAtom2_ (# of pairs) encountered since instantiation.
      int    maxNAtom2_;     
   
      /// The number of times this PairList has been built since instantiation.
      int    buildCounter_;
  
      /* 
      * Implementation Notes:
      *
      * The pair list is stored in Atom* pointer arrays atom1Ptrs_, atom2Ptrs_, 
      * and in an integer array first_. Each element of atom1Ptrs_ contains a 
      * pointer to the first atom in a pair (the primary Atom). Each element 
      * of atom2Ptrs_ contains a pointer to the second atom in a pair (the 
      * secondary Atom).  Secondary atoms that are neighbors of the same 
      * primary atom are listed consecutively.  Element first_[i] contains the 
      * array index of the first element in atom2Ptrs_ that contains a neighbor 
      * of primary atom atom1Ptrs_[i]. Pointers to neighbors of primary atom
      * *atom1Ptrs_[i] are thus stored in elements first_[i] <= j < first_[i+1]
      * of atom2Ptrs_.  Each pair is included only once, by including only 
      * pairs in which atom2Ptrs_[j]->id() > atom1Ptrs_[i]->id(), so
      * the Id of the secondary atom in a pair is always greater than that of
      * the primary atom.  The total number of primary atoms in atom1Ptrs_ is 
      * nAtom1_. The total number of secondary atoms in Array atom2Ptrs_ (or 
      * the total number of distinct pairs) is nAtom2_. The total number of
      * Atoms in the System in nAtom_.
      *
      * The only way legal way to loop over all atom pairs, using the public 
      * interface of a PairList, is to use a PairListIterator. See the 
      * documentation of PairListIterator for a discussion of its usage. The 
      * resulting loop over pairs is equivalent to the following double loop 
      * over private member variables of PairList:
      * 
      *    int i, j;
      *    Atom* atom1Ptr;
      *    Atom* atom2Ptr;
      *
      *    \\ Loop over primary Atoms
      *    for (i = 0; i < nAtom1_; ++i) {
      *       atom1Ptr = atom1Ptrs_[i];
      *
      *       // Loop over secondary atoms
      *       for (j = first_[i]; j < first_[i+1]; ++j) {
      *           atom2Ptr = atom2Ptrs_[j];
      *
      *           // ... Do something with Atoms *atom1Ptr and *atom2Ptr.
      *
      *       }
      *    }
      *
      * Note that DArray<int> first_ contains nAtom1_ + 1 elements. The first
      * element, element first_[0] is equal to 0, and the last element,
      * first_[nAtom1_], is always equal to nAtom2_, the total number of 
      * pairs.
      *
      * The number of primary Atoms nAtom1_ is generally less than nAtom_, the 
      * total number of atoms, because some Atoms have no neighbors with integer
      * Ids higher than their own. This is always true, for example, for the
      * Atom with the highest Id. Pointers to such atoms, which do not appear
      * in the first nAtom1_ elements of atom1Ptrs_, are instead stored  in the
      * last elements of the atom1Ptrs_ array , in elements tList1 + 1, .... ,
      * atomCapacity_ - 1. The previous position Vectors for these Atoms are 
      * stored in the corresponding elements of the oldPositions_ array. 
      * These pointers and previous positions retained because they are used 
      * in the isCurrent() method to calculate how far these Atoms have moved 
      * since the last time the PairList was built.
      */

   // friends:
   
      // Unit test classes
      friend class ::PairListTest;
   
   }; // end class PairList


   /*
   * Add an Atom to the CellList.
   */ 
   inline void PairList::addAtom(Atom &atom)
   { cellList_.addAtom(atom); }

   /*
   * Get the current number of atoms in the pairlist.
   */ 
   inline int PairList::nAtom() const
   { return nAtom_; }

   /*
   * Get the current number of pairs.
   */ 
   inline int PairList::nPair() const
   { return nAtom2_; }

   /*
   * Get the maximum value of aAtom() since instantiation.
   */ 
   inline int PairList::maxNAtom() const
   { return maxNAtom_; }

   /*
   * Get the maximum value of nPair() since instantiation.
   */ 
   inline int PairList::maxNPair() const
   { return maxNAtom2_; }

   /*
   * Get the number of times this PairList has been built.
   */ 
   inline int PairList::buildCounter() const
   { return buildCounter_; }

   /*
   * Has memory been allocated for this PairList?
   */ 
   inline bool PairList::isAllocated() const
   { return cellList_.isAllocated(); }

} 
#endif
