#ifndef DDMD_PAIR_LIST_H
#define DDMD_PAIR_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CellList.h"
#include <util/containers/GArray.h>
#include <util/misc/Setable.h>
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   class Atom;
   class PairIterator;
   
   /**
   * A Verlet nonbonded pair list.
   *
   * A PairList (or Verlet list) is a list of neighboring pairs of 
   * Atoms that are separated by a distance less than a specified 
   * cutoff distance. The PairList includes all pairs involving 
   * local atoms, i.e., all local-local or local-ghost pairs, with
   * each such pair included only once.
   *
   * To build or rebuild a PairList, one must first build an 
   * associated CellList, and then call PairList::build() to build 
   * the actual PairList. 
   *
   * A PairIterator object must be used to iterate over all of the 
   * pairs in in completed PairList (see documentation of PairIterator 
   * for usage). The first atom in each pair is always a local atom
   * and the second may be either a local or ghost atom.
   *
   * \ingroup DdMd_Neighbor_Module
   */
   class PairList 
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
      */
      virtual ~PairList();
   
      /// \name Mutators
      //@{

      /**
      * Set the pair list cutoff distance.
      *
      * May be called more than once, to reset the cutoff.
      *
      * \param cutoff pair list cutoff = potential cutoff + skin
      */
      void setCutoff(double cutoff);

      /**
      * Reserve space for atomCapacity atoms.
      *
      * \param atomCapacity maximum number of primary atoms
      */
      void reserveAtoms(int atomCapacity);

      /**
      * Reserve space for pairCapacity pairs.
      *
      * \param pairCapacity maximum number of pairs
      */
      void reservePairs(int pairCapacity);

      /**
      * Reset this to empty state.
      */  
      void clear();

      /**
      * Count and return the number of pairs.
      *
      * This function can be use to count the number of pairs before 
      * a simulation is begun, so that invoking code can choose an 
      * appropriate initial value for pairCapacity. The implementation
      * is basically a dry run of the algorithm used in the build()
      * function.
      *
      * \param cellList  CellList object that was just built
      * \param reverseUpdateFlag  True if reverse communication enabled
      * \return number of pairs within the pair list cutoff
      */
      int countPairs(CellList& cellList, bool reverseUpdateFlag = false);

      /**
      * Use a CellList to build a new PairList.
      *
      * \param cellList      a CellList object that was just built.
      * \param reverseUpdateFlag is reverse communication enabled?
      */
      void build(CellList& cellList, bool reverseUpdateFlag = false);

      //@}
      /// \name Accessors (miscellaneous)
      //@{ 
 
      /**
      * Initialize a PairIterator.
      *
      * \param iterator a PairList, initialized on output
      */
      void begin(PairIterator &iterator) const;
 
      /**
      * Get the number of primary atoms in the PairList.
      */
      int nAtom() const;

      /**
      * Get the number of pairs in the PairList.
      */
      int nPair() const;

      /**
      * Get the maximum number of pairs. 
      */
      int pairCapacity() const;

      /**
      * Get the maximum number of primary atoms.
      */
      int atomCapacity() const;

      //@}
      /// \name Statistics
      //@{
 
      /**
      * Compute statistics (reduce from all processors).
      * 
      * Call on all processors.
      */
      #ifdef UTIL_MPI
      virtual void computeStatistics(MPI::Intracomm& communicator);
      #else
      virtual void computeStatistics();
      #endif

      /**
      * Clear statistical accumulators (call on all processors).
      */
      void clearStatistics();

      /**
      * Output statistics.
      *
      * Call on master, after calling computeStatistics on all procs.
      *
      * \param out   output stream
      */
      void outputStatistics(std::ostream& out);

      /**
      * Get the maximum number of primary atoms encountered thus far.
      *
      * Call only on master.
      */
      int maxNAtom() const;

      /**
      * Get the maximum number of pairs encountered thus far.
      *
      * Call only on master.
      */
      int maxNPair() const;

      /**
      * Return number of times the PairList has been built thus far.
      */
      int buildCounter() const;

      //@}

   private:
  
      /// Array of pointers to 1st (primary) atom in each pair.
      GArray<Atom*>  atom1Ptrs_;  

      /// Array of pointers to neighbor (secondary) atom in each pair.
      GArray<Atom*>  atom2Ptrs_;  

      /// Array of indices in atom2Ptrs_ of first neighbor of an Atom.
      GArray<int>  first_; 

      /// Pair list cutoff radius (pair potential cutoff + skin_).
      double cutoff_;
   
      /// Maximum number of primary atoms on this proc since cleared.
      int  maxNAtomLocal_;     
   
      /// Maximum # of pairs on this proc since stats last cleared.
      int  maxNPairLocal_;     
   
      /// The number of times list has been built since stats cleared.
      int  buildCounter_;

      /// Max. number of primary atoms on all procs (defined only on master).
      Setable<int>  maxNAtom_;     
   
      /// Max. number of pairs on all procs (defined only on master).
      Setable<int>  maxNPair_;     
   
      /* 
      * Implementation Notes:
      *
      * The pair list is stored in Atom* pointer arrays atom1Ptrs_ and
      * atom2Ptrs_, and in an integer array first_. Each element of 
      * atom1Ptrs_ contains a pointer to the first atom in a pair (the 
      * primary Atom). Each element of atom2Ptrs_ contains a pointer 
      * to the second atom in a pair (the secondary Atom).  Secondary 
      * atoms that are neighbors of the same primary atom are listed 
      * consecutively.  Element first_[i] contains the array index of 
      * the first element in atom2Ptrs_ that contains a neighbor 
      * of primary atom atom1Ptrs_[i]. Pointers to neighbors of 
      * primary atom *atom1Ptrs_[i] are thus stored in elements 
      * first_[i] <= j < first_[i+1] of atom2Ptrs_.  Each pair is 
      * included only once. Primary atoms are all local atoms (i.e.,
      * are owned by this processor), while secondary atoms can be
      * local or ghost atoms.
      *
      * The only way legal way to loop over all atom pairs, using the 
      * public interface of a PairList, is to use a PairListIterator. 
      * See the documentation of PairListIterator for a discussion of 
      * its usage. The resulting loop over pairs is equivalent to the 
      * following double loop over private member variables of 
      * PairList:
      * 
      *    int i, j;
      *    Atom* atom1Ptr;
      *    Atom* atom2Ptr;
      *
      *    \\ Loop over primary Atoms
      *    for (i = 0; i < atom1Ptrs_.size(); ++i) {
      *       atom1Ptr = atom1Ptrs_[i];
      *
      *       // Loop over secondary atoms
      *       for (j = first_[i]; j < first_[i+1]; ++j) {
      *           atom2Ptr = atom2Ptrs_[j];
      *
      *           // Do something with Atoms *atom1Ptr and *atom2Ptr.
      *
      *       }
      *    }
      *
      * Note that GArray<int> first_ contains one more element than 
      * atom1Ptrs_. The element first_[0] is equal to 0, and the last 
      * element is always equal to the total number of pairs.
      */

   }; 

   /*
   * Get the current number of primary atoms in the pairlist.
   */ 
   inline int PairList::nAtom() const
   {  return atom1Ptrs_.size(); }

   /*
   * Get the current number of pairs.
   */ 
   inline int PairList::nPair() const
   {  return atom2Ptrs_.size(); }

   /*
   * Get the maximum number of pairs. 
   */
   inline int PairList::pairCapacity() const
   {  return atom2Ptrs_.capacity(); }

   /*
   * Get the maximum number of primary atoms.
   */
   inline int PairList::atomCapacity() const
   {  return atom1Ptrs_.capacity(); }

   /*
   * Get the maximum value of aAtom() since instantiation.
   */ 
   inline int PairList::maxNAtom() const
   { return maxNAtom_.value(); }

   /*
   * Get the maximum value of nPair() since instantiation.
   */ 
   inline int PairList::maxNPair() const
   { return maxNPair_.value(); }

   /*
   * Get the number of times this PairList has been built.
   */ 
   inline int PairList::buildCounter() const
   { return buildCounter_; }

} 
#endif
