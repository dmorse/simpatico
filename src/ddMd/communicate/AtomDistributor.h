#ifndef DDMD_ATOM_DISTRIBUTOR_H
#define DDMD_ATOM_DISTRIBUTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/chemistry/AtomArray.h>        // member
#include <ddMd/chemistry/Atom.h>             // member
#include <util/containers/ArrayStack.h>      // member
#include <util/containers/DMatrix.h>         // member
#include <util/containers/DArray.h>          // member
#include <util/boundary/Boundary.h>          // typedef

namespace DdMd
{

   class Domain;
   class Buffer;
   class AtomStorage;

   using namespace Util;

   /**
   * Class for distributing Atoms among processors.
   *
   * A AtomDistributor is used to distribute items among processors during
   * startup, when the master process must read a configuration file.
   *
   * Usage:
   * \code
   *
   *    AtomDistributor distributor;
   *
   *    Domain domain;
   *    Boundary boundary;
   *    AtomStorage storage;
   *    Buffer buffer;
   *    associate(domain, boundary, atomStorage, buffer);
   *
   *    if (rank = 0) {  // If master processorz
   *
   *       Atom* ptr;
   *       std::ifstream file
   *
   *       distributor.setup();
   *
   *       // Read from file
   *       for (i = 0; i < nAtom; ++i) {
   *
   *           ptr = distributor.newAtomPtr();
   *
   *           // Read properties of Atom *ptr from file
   *           file >> ptr->position();
   *           // ...
   *
   *           // Cache active atom *ptr for sending.
   *           distributor.addAtom(storage);
   *       }
   *
   *       // Send all remaining.
   *       distributor.send();
   *
   *    } else { // If not master processor
   *
   *       distributor.receive(storage);
   *
   *    }
   *
   * \endcode
   * The addAtom(storage) method can send items if required by
   * memory limits, and the send() method then sends all remaining
   * atoms.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class AtomDistributor : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      AtomDistributor();

      /**
      * Destructor.
      */
      ~AtomDistributor();

      /**
      * Set pointers to associated objects.
      *
      * This method must be called on all nodes, before any other.
      *
      * \param boundary  Boundary object (periodic boundary conditions)
      * \param domain    Domain object (processor grid)
      * \param storage   AtomStorage object (Atom container)
      * \param buffer    Buffer object (for communication)
      */
      void associate(Domain& domain, Boundary& boundary,
                     AtomStorage& storage, Buffer& buffer);

      /**
      * Set capacity of the send cache on the master node.
      *
      * This method may be called on master processor before the first
      * invokation of setup. Calling it on other processors has no 
      * effect, but does no harm. 
      *
      * Setting cacheCapacity < 0 enables computation of a default value 
      * large enough to accomodate full send buffers for all processors 
      * simultaneously.
      *
      * \param cacheCapacity max number of atoms cached for sending
      */
      void setCapacity(int cacheCapacity);

      /**
      * Read cacheCapacity.
      *
      * \param in input stream from which parameter is read.
      */
      virtual void readParameters(std::istream& in);

      #ifdef UTIL_MPI
      /**
      * Initialization before the loop over atoms on master processor.
      *
      * This method may only be called on the master processor,
      * just before entering the loop to read atoms from file. It
      * must be called after associate(), and after the associated
      * Domain and Buffer objects are initialized.
      */
      void setup();
      #endif

      /**
      * Returns pointer an address available for a new Atom.
      *
      * This method should be called only by the master processor. It
      * returns the address within the internal cache for a new Atom
      * object. Every call to newAtomPtr() must be followed by a 
      * matching call to addAtom().
      *
      * \return address for a new Atom.
      */
      Atom* newAtomPtr();

      /**
      * Process the active atom for sending.
      *
      * This method may be called only by the master processor, after
      * a matching call to newAtomPtr(). It identifies which processor
      * owns the active atom, i.e., the atom returned by the most recent
      * call to newAtomPtr. After shifting the position of this atom to
      * lie within the primary cell, it identifies the rank of the
      * processor that owns this atom, based on its position.  If this
      * atom is owned by the master (rank == 0), the atom is added to
      * the AtomStorage on the master node. Otherwise, the atom is added
      * to a cache of stored atoms, and a pointer to this atom is added
      * to a sendList for the relevant processor.
      *
      * Atoms are actually transmitted to other processors as needed:
      * Whenever addition of an atom to the cache would exceed either
      * the cache capacity or the capacity of the sendlist for the owner
      * processor, this method sends a buffer to the processor with the
      * the largest sendList, before caching the current atom and marking
      * it for later sending.
      *
      * The coordinates of the active atom must be expressed in 
      * scaled / generalized [0, 1] coordinates when it is added.
      *
      * \return rank of processor that owns the active atom.
      */
      int addAtom();

      #ifdef UTIL_MPI
      /**
      * Send all atoms that have not be sent previously.
      *
      * This method should be called only by the master processor, after
      * completion of a loop over all atoms to be sent.
      */
      void send();

      /**
      * Receive all atoms sent by master processor.
      *
      * This should be called by all processes except the master. 
      *
      * Upon return, all processors should have the correct atoms, 
      * with scaled / generalized atom coordinates set correctly.
      */
      void receive();
      #endif

      /**
      * Validate distribution of atoms after completion.
      *
      * This method checks call AtomStorage::isValid(), checks that
      * the total number of atoms on all processors is equal to nAtom,
      * and checks that every atom is in the correct processor domain.
      *
      * On the master processor
      * \return total number of atoms distributed on master, 0 on slaves.
      */
      int validate();

   private:

      /**
      * Array that holds cached Atom objects to be sent to other processors.
      * Allocated only on the master processor.
      */
      AtomArray cache_;

      /**
      * Stack of pointers to unused elements of cache_ array.
      * Allocated only on the master processor.
      */
      ArrayStack<Atom> reservoir_;

      #ifdef UTIL_MPI
      /**
      * Matrix of ptrs to elements of cache_ array. Each row contains
      * pointers to atoms to be sent to one processor.
      * Allocated only on the master processor.
      */
      DMatrix<Atom*> sendArrays_;

      /**
      * Array of sendArrays_ row sizes.
      * Allocated only on the master processor.
      */
      DArray<int> sendSizes_;
      #endif

      /**
      * Pointer to associated Boundary object (all domain processors).
      */
      Boundary* boundaryPtr_;

      /**
      * Pointer to associated Domain object (all domain processors).
      */
      Domain* domainPtr_;

      /**
      * Pointer to associated Domain object (all domain processors).
      */
      AtomStorage* storagePtr_;

      #ifdef UTIL_MPI
      /**
      * Pointer to associated Buffer object (all domain processors).
      */
      Buffer* bufferPtr_;
      #endif

      /**
      * Pointer to space for a new local Atom. Null when inactive.
      */
      Atom* newPtr_;

      /**
      * Maximum number of items that can be cached on master.
      */
      int cacheCapacity_;

      /**
      * Maximum allowed number of items per transmission to one processor.
      */
      int sendCapacity_;

      /**
      * Rank of processor with the maximum current buffer size.
      * Used only on the master processor.
      */
      int rankMaxSendSize_;

      /**
      * Total number of items cached for sending in addAtom()
      * Used only on the master processor.
      */
      int nCachedTotal_;

      /**
      * Total number of items actually sent.
      * Used only on the master processor.
      */
      int nSentTotal_;

      /**
      * Has memory been allocated?
      *
      * Initialized to false, set true at end of allocate() method.
      */
      bool isAllocated_;

      /**
      * Allocate and initialize memory on master processor (private).
      *
      * Only called on the master processor. This function allocates 
      * and initializes all required memory on the master.  It is 
      * called by setup on the first invokation of setup. It assigns
      * isAllocated_ = true before returning. 
      */
      void allocate();

   };

}
#endif
