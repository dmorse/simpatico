#ifndef DDMD_ATOM_COLLECTOR_H
#define DDMD_ATOM_COLLECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/chemistry/AtomArray.h>        // member
#include <ddMd/storage/AtomIterator.h>       // member

namespace DdMd
{

   class AtomStorage;
   class Domain;
   class Buffer;

   using namespace Util;

   /**
   * Class for collecting Atoms from processors to master processor.
   *
   * An AtomCollector collects atoms from all processors to the master
   * in order to output a configuration file. 
   *
   * Usage:
   * \code
   * 
   *    AtomStorage storage;
   *    Domain domain;
   *    Buffer buffer;
   *    AtomCollector collector;
   *
   *    // Initialization
   *    collector.associate(domain, storage, buffer);
   *
   *    // Communication
   *    if (domain.gridRank() == 0) {  // if master processor
   *       collector.setup();
   *       Atom* atomPtr = collector.nextPtr();
   *       while (atomPtr) {
   *          // Write *atomPtr to file; 
   *          atomPtr = collector.nextPtr();
   *       }
   *    } else { // if not master
   *       collector.send();
   *    }
   *
   * \endcode
   *
   * \ingroup DdMd_Communicate_Module
   */
   class AtomCollector : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      AtomCollector();

      /**
      * Destructor.
      */
      ~AtomCollector();

      /**
      * Initialize pointers to associated objects.
      *
      * Call on all processors, only once, before first call to setup.
      */
      void associate(Domain& domain, AtomStorage& storage, Buffer& buffer);

      /**
      * Set cache capacity on master processor.
      *
      * Call oly once on master processors, before first call to setup.
      * Calling on other processors has no effect, but does no harm.
      *
      * \param recvArrayCapacity size of array used to receive atoms.
      */
      void setCapacity(int recvArrayCapacity);

      /**
      * Setup master processor for receiving.
      *
      * Call only on the master processor, just before each receive loop.
      */
      void setup();

      /**
      * Return a pointer to the next available atom, or null. 
      *
      * Call only on the master processor, within loop over atoms.
      *
      * \return address of next atom, or null if no more are available.
      */
      Atom* nextPtr();
     
      /**
      * Send all atoms to the master.
      *
      * Call on all processors except the master.
      */
      void send();

   private:

      /// Temporary array of atoms, allocated only on master.
      AtomArray recvArray_;

      /// Iterator for atoms in a AtomStorage (on all domain nodes).
      AtomIterator iterator_;

      /// Pointer to associated Domain object (on all domain nodes).
      Domain* domainPtr_;

      /// Pointer to associated AtomStorage object (on all domain nodes).
      AtomStorage* storagePtr_;

      /// Pointer to associated Buffer object (on all domain nodes).
      Buffer* bufferPtr_;

      /// Rank of processor from which atoms are being received (on master).
      int source_;

      /// Number of items in receive buffer (on master).
      int recvBufferSize_;

      /// Capacity of recvArray_ (allocated only on master, known on all).
      int recvArrayCapacity_;

      /// Number of items in recvArray_ (on master).
      int recvArraySize_;

      /// Index of current item in recvArray_ (on master).
      int recvArrayId_;

      /// Have all atoms been processed from current source? (all domain nodes).
      bool isComplete_;

   };

}
#endif
