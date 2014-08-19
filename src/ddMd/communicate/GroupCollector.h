#ifndef DDMD_GROUP_COLLECTOR_H
#define DDMD_GROUP_COLLECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/storage/GroupIterator.h>      // member
#include <ddMd/chemistry/Group.h>            // member
#include <util/containers/DArray.h>          // member

namespace DdMd
{

   class Domain;
   class Buffer;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Class for collecting Groups from processors to master processor.
   *
   * An GroupCollector collects Group objects from all processors to
   * the master processor in order to output a configuration file. 
   *
   * Usage:
   * \code
   * 
   *    GroupStorage<N> storage;
   *    Domain domain;
   *    Buffer buffer;
   *    GroupCollector collector;
   *
   *    // Initialization
   *    collector.associate(domain, storage, buffer);
   *
   *    // Communication
   *    if (domain.gridRank() == 0) {  // if master processor
   *       collector.allocate(100);
   *       collector.setup();
   *       Group<N>* groupPtr = collector.nextPtr();
   *       while (groupPtr) {
   *          // Write *groupPtr to file; 
   *          groupPtr = collector.nextPtr();
   *       }
   *    } else { // if not master
   *       collector.send();
   *    }
   *
   * \endcode
   *
   * \ingroup DdMd_Communicate_Module
   */
   template <int N>
   class GroupCollector : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      GroupCollector();

      /**
      * Destructor.
      */
      ~GroupCollector();

      /**
      * Initialize pointers to associated objects.
      *
      * Call on all processors, only once.
      */
      void associate(Domain& domain, GroupStorage<N>& storage, Buffer& buffer);

      /**
      * Set size of cache for receiving groups on the master.
      *
      * The cache is actually allocated only on  master processor, within
      * the first call to the setup function. This function must be called 
      * once on the master. Calling it on other processors sets an unused
      * variable, with no other effect, and thus does no harm. 
      *
      * \param recvArrayCapacity capacity of recvArray cache on master.
      */
      void setCapacity(int recvArrayCapacity);

      /**
      * Setup master processor for receiving.
      *
      * Call only on the master processor, just before receive loop.
      */
      void setup();

      /**
      * Return a pointer to the next available atom, or null. 
      *
      * Call only on the master processor, within loop over groups. This
      * function return the address of the next available Group, or null
      * if no more are available. Groups are returned in order of 
      * processor rank, starting with those on the master processor, 
      * then processor 1, etc. 
      *
      * The address returned is an address within an internal cache, 
      * which is then freed for reuse. All Group data must thus be 
      * copied to a permanent location before calling nextPtr() again.
      *
      * \return address of next group, or null if no more are available.
      */
      Group<N>* nextPtr();
     
      /**
      * Send all groups on this processor to the master processor.
      *
      * Call on all processors except the master (rank = 0) processor.
      */
      void send();

   private:

      /// Temporary array of groups. Allocated only on master.
      DArray< Group<N> > recvArray_;

      /// Iterator for groups in a GroupStorage<N> (on all procs).
      GroupIterator<N> iterator_;

      /// Pointer to associated Domain object (on master).
      Domain* domainPtr_;

      /// Pointer to associated Domain object (on master).
      GroupStorage<N>* storagePtr_;

      /// Pointer to associated Buffer object (on master).
      Buffer* bufferPtr_;

      /// Rank of processor from which groups are being received (on master).
      int source_;

      /// Capacity of recvArray cache (allocated on master).
      int recvArrayCapacity_;

      /// Number of unread groups in MPI receive buffer (on master).
      int recvBufferSize_;

      /// Number of items in recvArray_ (on master).
      int recvArraySize_;

      /// Index of current item in recvArray_ (on master).
      int recvArrayId_;

      /// Have all groups been processed from current source? (all).
      bool isComplete_;

   };

}
#endif
