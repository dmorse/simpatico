#ifndef DDMD_GROUP_DISTRIBUTOR_H
#define DDMD_GROUP_DISTRIBUTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/chemistry/Group.h>            // member
#include <ddMd/communicate/Buffer.h>         // member data type
#include <util/containers/DArray.h>          // member

namespace DdMd
{

   class Domain;
   class AtomStorage;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Class template for distributing Group<N> objects among processors.
   *
   * A GroupDistributor is used to distribute items among processors during
   * startup, when the master process must read a configuration file. Each
   * ConfigIo has a GroupDistributor for each group type.
   *
   * Precondition:
   *
   * The usage pattern described below must be invoked after all local atoms
   * are distributed, but when there are no ghosts.  
   *
   * Usage:
   * \code
   * 
   *    GroupDistributor  distributor;
   *    GroupStorage<N>   groupStorage;
   *    Domain&           domain;
   *    AtomStorage       atomStorage, 
   *    Buffer            buffer;
   *    distributor.associate(Domain, atomStorage, GroupStorage, buffer);
   *
   *    DArray<int>       atomOwners;
   *    Group<N>*         ptr;
   *
   *    if (rank = 0) {  // If master processor
   *
   *       distributor.setup();
   *
   *       // Read from file
   *       std::ifstream file
   *       for (i = 0; i < nGroup; ++i) {
   *
   *           ptr = distributor.newPtr();
   *  
   *           // Read properties of Group<N> *ptr from file
   *           file >> ptr->position();
   *           // ...
   *
   *           // Cache active Group<N> for sending.
   *           distributor.add();
   *       }
   *
   *       // Send all remaining.
   *       distributor.send();
   *
   *    } else { // If not master processor
   *
   *       distributor.receive();
   *
   *    }
   * \endcode
   * The add(groupStorage) method can send items if required by memory
   * limits, and the send() method then sends all remaining groups.
   *
   * Posconditions:
   * 
   *   - At the end end of the above pattern, each processor has a copy 
   *     of every Group that contains one or more Atoms that it owns. 
   *
   *   - All pointers to local atoms are set, and pointer to missing atoms
   *     are set to null. Missing atoms will later be assigned to ghosts.
   *
   * \ingroup DdMd_Communicate_Module
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
      * Create required associations with related objects.
      *
      * \param domain        associated Domain object defines the processor grid
      * \param atomStorage   associated AtomStorage manages memory for atoms
      * \param groupStorage  associated GroupStorage manages memory for groups
      * \param buffer        associated buffer provides memory for communication
      */
      void associate(Domain& domain, 
                     AtomStorage& atomStorage, 
                     GroupStorage<N>& groupStorage, 
                     Buffer& buffer);

      /**
      * Initialize Buffer for sending.
      */
      void setup();

      /**
      * Set cacheCapacity, allocate memory and initialize object.
      *
      * This method sets cacheCapacity from file, and then goes through
      * the same initialization steps as readParam.
      *
      * \param cacheCapacity max number of groups cached for sending
      */
      void setCapacity(int cacheCapacity);

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
      virtual void readParameters(std::istream& in);

      /**
      * Returns pointer an address available for a new Group<N>.
      *
      * This method should be called only by the master processor. It
      * returns the address within the internal cache for a new Group<N>
      * object. Each call to newPtr() must be followed by a matching
      * call to add(). 
      *
      * \return address for a new Group<N> object.
      */
      Group<N>* newPtr(); 
     
      /**
      * Add a group to the cache for sending, send if necessary.
      *
      * If the addition of this atom to the cache would make the cache
      * full, this method broadcasts the cache to all processors.
      */
      void add();

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
      void receive();

   protected:

      /**
      * Allocate memory and initialize object.
      *
      * This method allocates all required memory and leaves the object ready
      * for use. It is called by setParam and readParam.
      */
      void allocate();

   private:

      /// Array of cached Group objects to be broadcast to all other procs.
      /// Allocated only on the master processor.
      DArray< Group<N> > cache_;

      /// Pointer to space in cache for a new local Group<N>. Null when inactive.
      Group<N>* newPtr_;
      
      /// Pointer to associated Domain object.
      Domain* domainPtr_;

      /// Pointer to associated AtomStorage object.
      AtomStorage* atomStoragePtr_;

      /// Pointer to associated GroupStorage<N> object.
      GroupStorage<N>* groupStoragePtr_;

      /// Pointer to associated Buffer object.
      Buffer* bufferPtr_;

      /// Type of object to send.
      enum Buffer::BlockDataType sendType_;

      /// Total number of atoms in groups recieved (defined on all).
      int nAtomRecv_;

      /// Total number of groups sent (defined only on master).
      int nSentTotal_;

      /// Allocated capacity of cache_ (allocated only on master).
      int cacheCapacity_;

      /// Current size of cache_ (defined only on master).
      int cacheSize_;

      /**
      * Validate groups after receipt.
      */
      void validate();

   };

}
#endif
