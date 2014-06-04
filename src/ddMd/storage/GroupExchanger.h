#ifndef DDMD_GROUP_EXCHANGER_H
#define DDMD_GROUP_EXCHANGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>             
#include <util/space/Dimension.h> 

namespace Util {
   class IntVector;
   template <typename T, int M, int N> class FMatrix;
   template <typename T> class GPArray; 
}

namespace DdMd
{

   class Atom;
   class AtomStorage;
   class Buffer;
   using namespace Util;

   /**
   * Interface for a GroupStorage<N> for use in Exchanger.
   *
   * \ingroup DdMd_Storage_Module
   */
   class GroupExchanger
   {

   public:

      /**
      * Default constructor.
      */
      GroupExchanger();

      /**
      * Destructor.
      */
      virtual ~GroupExchanger();

      /**
      * Find and mark groups that span boundaries.
      *
      * \param bound  boundaries of domain for this processor
      * \param inner  inner slab boundaries (extended domain of neighbors)
      * \param outer  outer slab boundaries (extended domain of this processor)
      * \param gridFlags  element i is 0 iff gridDimension[i] == 1, 1 otherwise
      */
      virtual
      void markSpanningGroups(FMatrix<double, Dimension, 2>& bound, 
                              FMatrix<double, Dimension, 2>& inner, 
                              FMatrix<double, Dimension, 2>& outer, 
                              IntVector& gridFlags) = 0;
   
      #ifdef UTIL_MPI
      /**
      * Pack groups for exchange.
      *
      * Usage: This is called after atom exchange plans are set, 
      * but before atoms marked for exchange in direction i, j
      * have been removed from the atomStorage.
      *
      * \param i  index of Cartesian communication direction
      * \param j  index for sign of direction
      * \param buffer  Buffer object into which groups are packed
      */
      virtual void pack(int i, int j, Buffer& buffer) = 0;
   
      /**
      * Unpack groups from buffer and find available associated atoms.
      *
      * This function should unpack groups, add new ones to a GroupStorage,
      * set pointers to all Group atoms that are in the AtomStorage, and
      * nullify pointers to atoms that are not present.
      *
      * \param buffer  Buffer object from which groups are unpacked
      * \param atomStorage  AtomStorage used to find atoms pointers
      */
      virtual void unpack(Buffer& buffer, AtomStorage& atomStorage) = 0;
      #endif // ifdef UTIL_MPI
   
      /**
      * Set ghost communication plans for all atoms in incomplete groups.
      *
      * Usage: This is called after exchanging all atoms and groups between 
      * processors, but before exchanging ghosts. At this point, atom
      * ownership is finalized, but there are no ghosts.
      *
      * \param atomStorage AtomStorage object
      * \param sendArray   Matrix of arrays of pointers to ghosts to send
      * \param gridFlags   element i is 0 iff gridDimension[i] == 1, 1 otherwise
      */
      virtual
      void markGhosts(AtomStorage& atomStorage, 
                      FMatrix< GPArray<Atom>, Dimension, 2 >&  sendArray,
                      IntVector& gridFlags) = 0;
   
      /**
      * Find all ghost members of groups.
      *
      * Usage: This called after all ghosts have been exchanged.
      *
      * \param atomStorage AtomStorage object used to find atom pointers
      */
      virtual void findGhosts(AtomStorage& atomStorage) = 0;
   
      /**
      * Return true if the container is valid, or throw an Exception.
      *
      * Calls overloaded isValid() function, then checks consistency of atom 
      * pointers with those in an asociated AtomStorage. If hasGhosts is
      * false, the function requires that no group contain a pointer to a ghost 
      * atom. If hasGhosts is true, requires that every Group be complete.
      *
      * \param atomStorage  associated AtomStorage object
      * \param hasGhosts    true if the atomStorage has ghosts, false otherwise
      * \param communicator domain communicator 
      */
      #ifdef UTIL_MPI
      virtual 
      bool isValid(AtomStorage& atomStorage, MPI::Intracomm& communicator, 
                   bool hasGhosts) = 0;
      #else
      virtual 
      bool isValid(AtomStorage& atomStorage, bool hasGhosts) = 0;
      #endif

   }; // class GroupExchanger

} // namespace DdMd
#endif
