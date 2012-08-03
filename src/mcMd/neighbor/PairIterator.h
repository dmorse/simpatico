#ifndef MCMD_PAIR_ITERATOR_H
#define MCMD_PAIR_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PairList.h"
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class Atom;
   
   /**
   * Iterator for pairs in a PairList. 
   *
   * A loop over atom pairs, such as that required to calculated forces or
   * energies, can take the following form:
   * \code
   * 
   *    PairList     pairList;
   *    PairIterator iter;
   *
   *    // [build the PairList]
   *
   *    // Loop over pairs
   *    Atom* atom1Ptr;
   *    Atom* atom2Ptr;
   *    for (pairList.begin(iter); iterator.notEnd(); ++iter) {
   *       iterator.getPair(atom1Ptr, atom2Ptr);
   *
   *       // [Do something with atom1Ptr and atom2Ptr]
   *
   *       ++iter;
   *    }
   *
   * \endcode
   *
   * \ingroup McMd_Neighbor_Module
   */
   class PairIterator
   {
   
   public:

      /**
      * Default constructor.
      *
      * Creates an uninitialized iterator, which is not associated with a 
      * parent PairList.
      */
      PairIterator();
   
      /**
      * Initializing constructor.
      *
      * Creates an initialized iterator, for which the current pair is the
      * first pair in a parent PairList.
      *
      * \param pairList parent PairList object.
      */
      PairIterator(const PairList &pairList);
  
      // Use default C++ destructor.
 
      /** 
      * Increment to next pair.
      */
      PairIterator& operator++ ();
 
      /** 
      * Return true if at end of PairList.
      *
      * When isEnd() is true, the current pair is already one past past
      * the end of the pairlist, and is thus invalid.
      */
      bool isEnd() const;
   
      /** 
      * Return true if not at end of PairList.
      *
      * When notEnd() is false, the current pair is already one past past
      * the end of the pairlist, and is thus invalid.
      */
      bool notEnd() const;
   
      /**
      * Get pointers for current pair of Atoms. 
      *
      * \param atom1Ptr pointer to current atom 1.
      * \param atom2Ptr pointer to current atom 2.
      */
      void getPair(Atom* &atom1Ptr, Atom* &atom2Ptr) const;

   private:
 
      /// Array of const pointers to primary atom in each pair.
      Atom* const * atom1Ptrs_;  

      /// Array of const pointers to secondary atom in each pair.
      Atom* const * atom2Ptrs_;  

      /// Pointer to const index in atom2Ptrs_ of first neighbor of an Atom.
      const int *   first_; 

      /// Number of primary atoms in atom1Ptrs_.
      int    nAtom1_;      
  
      /// Number of secondary atoms in atom2Ptrs_.
      int    nAtom2_;      
  
      /// Current index of first atom in atom1Ptrs_.
      int    atom1Id_;

      /// Current index of second atom in atom2Ptrs_.
      int    atom2Id_;

   //friend:

      friend class PairList;

   }; // end class PairIterator


   /*
   * Default constructor.
   */
   inline PairIterator::PairIterator()
    : atom1Ptrs_(0),
      atom2Ptrs_(0),
      first_(0),
      nAtom1_(0),
      nAtom2_(0),
      atom1Id_(0),
      atom2Id_(0)
   {}

   /*
   * Constructor, initialized iterator.
   */
   inline PairIterator::PairIterator(const PairList &pairList)
    : atom1Ptrs_(0),
      atom2Ptrs_(0),
      first_(0),
      nAtom1_(0),
      nAtom2_(0),
      atom1Id_(0),
      atom2Id_(0)
   {  pairList.begin(*this); }

   /*
   * Get pointers to current pair of Atoms. 
   */
   inline void PairIterator::getPair(Atom* &atom1Ptr, Atom* &atom2Ptr) const
   {
      assert(atom1Id_ >=0);
      assert(atom1Id_ < nAtom1_);
      assert(atom2Id_ >=0);
      assert(atom2Id_ < nAtom2_);
      atom1Ptr = atom1Ptrs_[atom1Id_];
      atom2Ptr = atom2Ptrs_[atom2Id_];
   }
 
   /*
   * Increment current pair.
   */
   inline PairIterator& PairIterator::operator++ ()
   {
      assert(atom1Id_ >=0);
      assert(atom1Id_ < nAtom1_);
      assert(atom2Id_ >=0);
      assert(atom2Id_ < nAtom2_);
      ++atom2Id_;
      if (atom2Id_ == first_[atom1Id_+1]) {
         ++atom1Id_;
      }
      return *this;
   }

   /*
   * Return true if at end of the pair list.
   */
   inline bool PairIterator::isEnd() const
   { return (atom2Id_ == nAtom2_); }
   
   /*
   * Return true if not at end of pair list.
   */
   inline bool PairIterator::notEnd() const
   { return (atom2Id_ != nAtom2_); }
   
 
} 
#endif
