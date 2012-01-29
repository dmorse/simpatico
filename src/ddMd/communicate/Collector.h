#ifndef COLLECTOR_H
#define COLLECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   * A Collector is used to collect atoms from processors back to the
   * master in order to output a configuration file. 
   *
   * Usage:
   * \code
   * 
   *    AtomStorage storage;
   *    Domain      domain;
   *    Buffer      buffer;
   *    Collector   collector;
   *
   *    // If master processor
   *    if (domain.gridRank() == 0) {  
   *
   *       collector.initialize(storage, domain, buffer);
   *       Atom* atomPtr = collector.nextPtr();
   *       while (atomPtr) {
   *          // Write *atomPtr to file; 
   *          atomPtr = collector.nextPtr();
   *       }
   *
   *    } else { 
   *       // If not master processor
   *
   *       collector.send(storage, domain, buffer);
   *
   *    }
   *
   * \endcode
   */
   class Collector : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Collector();

      /**
      * Destructor.
      */
      ~Collector();

      /**
      * Initialize master processor for receiving.
      *
      * Call only on the master processor.
      */
      void initialize(AtomStorage& storage, Domain& domain, Buffer& buffer);

      /**
      * Return a pointer to the next available atom, or null. 
      *
      * Call only on the master processor.
      *
      * \return address of next atom, or null if no more are available.
      */
      Atom* nextPtr();
     
      /**
      * Send all atoms to the master.
      *
      * Call on all processors except the master.
      */
      void send(AtomStorage& storage, Domain& domain, Buffer& buffer);

   private:

      /// Temporary array of atoms, allocated only on master.
      AtomArray   recvArray_;

      /// Iterator for atoms in a AtomStorage (master and slaves).
      AtomIterator iterator_;

      /// Pointer to associated Domain object (on master).
      Domain*     domainPtr_;

      /// Pointer to associated Buffer object (on master).
      Buffer*     bufferPtr_;

      /// Rank of processor from which atoms are being received (on master).
      int         source_;

      /// Number of items in recvArray_ (on master).
      int         recvArraySize_;

      /// Index of current item in recvArray_ (on master).
      int         recvArrayId_;

      /// Have all atoms been processed from current source? (master and slaves).
      bool        isComplete_;

   };

}
#endif
