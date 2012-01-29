#ifndef GROUP_DISTRIBUTOR_H
#define GROUP_DISTRIBUTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/chemistry/Group.h>            // member
#include <util/containers/DArray.h>          // member
#include <util/containers/DMatrix.h>         // member
#include <util/containers/ArrayStack.h>      // member

namespace DdMd
{

   class Domain;
   class Buffer;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Class template for distributing Group<N> objects among processors.
   *
   * A GroupDistributor is used to distribute items among processors during
   * startup, when the master process must read a configuration file.  The
   * loops to distribute Groups (bonds, angles, or dihedrals) requires access
   * to an DArray<int> atomOwners[i] in which atomOwners[i] is the rank of 
   * the processor that owns atom i.
   *
   * Usage:
   * \code
   * 
   *    GroupDistributor  distributor;
   *    GroupStorage<N>   storage;
   *    DArray<int>       atomOwners;
   *    Group<N>*         ptr;
   *    std::ifstream file
   *
   *    if (rank = 0) {  // If master processor
   *
   *       distributor.initSendBuffer();
   *
   *       // Read from file
   *       for (i = 0; i < nGroup; ++i) {
   *
   *           ptr = distributor.newPtr();
   *  
   *           // Read properties of Group<N> *ptr from file
   *           file >> ptr->position();
   *           // ...
   *
   *           // Cache active Group<N> for sending.
   *           distributor.add(storage, atomOwners);
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
   * The add(storage, atomOwners) method can send items if required by
   * memory limits, and the send() method then sends all remaining
   * groups.
   *
   */
   template <int N>
   class GroupDistributor : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      GroupDistributor();

      /**
      * Destructor.
      */
      ~GroupDistributor();

      /**
      * Set pointers to boundary, domain, and buffer.
      *
      * This method must be called on all nodes, before any other.
      * 
      * \param domain        Domain object (processor grid)
      * \param buffer        Buffer used for communication
      */
      void associate(Domain& domain, Buffer& buffer);

      /**
      * Set cacheCapacity, allocate memory and initialize object.
      *
      * This method sets cacheCapacity from file, and then goes through
      * the same initialization steps as readParam.
      *
      * \param cacheCapacity max number of groups cached for sending
      */
      void setParam(int cacheCapacity = -1);

      /**
      * Read cacheCapacity, allocate memory and initialize object.
      *
      * This method and setParam both allocate all required memory and 
      * initialize the GroupDistributor. Memory for a temporary read 
      * cache is alocated only on the master.
      * 
      * Preconditions: The initialize method must have been called, the 
      * Buffer must be allocated, and the Domain must be initialized.
      *
      * \param in input stream from which parameter is read.
      */
      virtual void readParam(std::istream& in);

      #ifdef UTIL_MPI

      /**
      * Initialize buffer before the loop over groups.
      *
      * This method should be called only by the master processor,
      * just before entering the loop to read groups from file. 
      */
      void initSendBuffer();

      #endif

      /**
      * Returns pointer an address available for a new Group<N>.
      *
      * This method should be called only by the master processor. It
      * returns the address within the internal cache for a new Group<N>
      * object.  Each call to newPtr() must be followed by a matching
      * call to add(). 
      *
      * \return address for a new Group<N> object.
      */
      Group<N>* newPtr(); 
     
      /**
      * Process the active group for sending.
      *
      * This method should be called only by the master processor, after
      * a matching call to newPtr(). It identifies and returns the
      * rank of the processor that owns the active group (i.e., the atom 
      * returned by the recent call to newPtr), based on its position. 
      * If this group is owned by the master (rank == 0), it adds it to
      * the GroupStorage<N> on the master node. If it is not owned by the
      * master, it caches the atom for sending to its owner. 
      *
      * If the addition of this atom to the cache would make the cache
      * full, this method first sends a buffer to the processor with
      * the largest sendList.
      *
      * \param  storage GroupStorage<N> object on master processor.
      * \return rank of processor that owns the active atom.
      */
      void add(GroupStorage<N>& storage, DArray<int> atomOwners);

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
      */ 
      void receive(GroupStorage<N>& storage);

   protected:

      /// Maximum number of items that can be cached on master.
      int         cacheCapacity_;

      /**
      * Allocate memory and initialize object.
      *
      * This method allocates all required memory and leaves the object ready
      * for use. It is called by setParam and readParam.
      */
      void allocate();

   private:

      /// Array that holds cached Group objects to be sent to other processors.
      /// Allocated only on the master processor.
      DArray< Group<N> > cache_;

      /// Stack of pointers to unused elements of cache_ array.
      /// Allocated only on the master processor.
      ArrayStack< Group<N> > reservoir_;
      
      /// Matrix of ptrs to elements of cache_ array. Each row contains
      /// pointers to Group objects to be sent to one processor.
      /// Allocated only on the master processor.
      DMatrix< Group<N>* >  sendArrays_;

      /// Array of sendArrays_ row sizes.
      /// Allocated only on the master processor.
      DArray<int> sendSizes_;

      /// Pointer to associated Domain object.
      Domain*     domainPtr_;

      /// Pointer to associated Buffer object.
      Buffer*     bufferPtr_;

      /// Pointer to space for a new local Group<N>. Null when inactive.
      Group<N>*   newPtr_;
      
      /// Maximum allowed number of items per transmission to one processor.
      int         sendCapacity_;

      /// Rank of processor with the maximum current buffer size.
      /// Used only on the master processor.
      int         rankMaxSendSize_;

      /// Total number of items cached for sending in add()
      /// Used only on the master processor.
      int          nCachedTotal_;

      /// Total number of items actually sent.
      /// Used only on the master processor.
      int          nSentTotal_;

   };

}
#endif
