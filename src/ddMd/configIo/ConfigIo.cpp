#ifndef DDMD_CONFIG_IO_CPP
#define DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

#include <ddMd/simulation/Simulation.h>                 
#include <ddMd/communicate/Domain.h>   
#include <ddMd/storage/AtomStorage.h>               
#include <ddMd/storage/BondStorage.h>               
#include <ddMd/communicate/GroupCollector_inc.h> 
#include <ddMd/communicate/GroupDistributor.cpp> 
#include <ddMd/communicate/Buffer.h> 
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Bond.h>
#include <ddMd/chemistry/MaskPolicy.h>
#include <util/space/Vector.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo()
    : domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0),
      bondStoragePtr_(0),
      atomCacheCapacity_(0),
      bondCacheCapacity_(0)
   {}

   /*
   * Read cache size and allocate memory.
   */
   void ConfigIo::associate(Domain& domain, Boundary& boundary,
                            AtomStorage& atomStorage,
                            BondStorage& bondStorage,
                            Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      atomStoragePtr_ = &atomStorage;
      bondStoragePtr_ = &bondStorage;
      atomDistributor_.associate(domain, boundary, atomStorage, buffer);
      bondDistributor_.associate(domain, atomStorage,
                                 bondStorage, buffer);
      atomCollector_.associate(domain, atomStorage, buffer);
      bondCollector_.associate(domain, bondStorage, buffer);
   }

   /*
   * Read cache size and allocate memory.
   */
   void ConfigIo::readParam(std::istream& in)
   {
      readBegin(in, "ConfigIo");
      read<int>(in, "atomCacheCapacity", atomCacheCapacity_);
      read<int>(in, "bondCacheCapacity", bondCacheCapacity_);
      readEnd(in);
      atomDistributor_.setParam(atomCacheCapacity_);
      bondDistributor_.setParam(bondCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
   }

   /*
   * Set parameters and allocate memory.
   */
   void ConfigIo::initialize(int atomCacheCapacity, int bondCacheCapacity)
   {
      atomCacheCapacity_ = atomCacheCapacity;
      bondCacheCapacity_ = bondCacheCapacity;
      atomDistributor_.setParam(atomCacheCapacity_);
      bondDistributor_.setParam(bondCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
   }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int ConfigIo::readGroups(std::istream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupDistributor<N>& distributor) 
   {
      int nGroup;  // Total number of groups in file
      if (domain().isMaster()) {  
         file >> Label(sectionLabel);
         file >> Label(nGroupLabel) >> nGroup;
         Bond* groupPtr;
         int i, j, k;
         distributor.initSendBuffer();
         for (i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            file >> *groupPtr;
            for (j = 0; j < 2; ++j) {
               k = groupPtr->atomId(j);
            }
            distributor.add();
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         // Receive all groups into BondStorage
         distributor.receive();
      }
      return nGroup;
   }

   /*
   * Read a configuration file.
   */
   void ConfigIo::readConfig(std::istream& file, MaskPolicy maskPolicy)
   {

      // Read and broadcast boundary
      if (domain().isMaster()) {  
         file >> Label("BOUNDARY");
         file >> boundary();
      }
      #if UTIL_MPI
      bcast(domainPtr_->communicator(), boundary(), 0);
      #endif

      // Atoms 
      int nAtom;  // Total number of atoms in file
      if (domain().isMaster()) {  

         // Read and distribute atoms
         file >> Label("ATOMS");

         // Read number of atoms
         file >> Label("nAtom") >> nAtom;

         int totalAtomCapacity = atomStoragePtr_->totalAtomCapacity();

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().initSendBuffer();
         #endif

         // Read atoms
         Atom* atomPtr;
         int id;
         int typeId;
         int rank;
         for(int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            file >> id >> typeId;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
            file >> atomPtr->position();
            file >> atomPtr->velocity();

            // Add atom to list for sending.
            rank = atomDistributor().addAtom();

         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor
         atomDistributor().receive();
      }

      // Check that all atoms are accounted for after distribution.
      {
         int nAtomLocal = atomStorage().nAtom();
         int nAtomAll;
         #ifdef UTIL_MPI
         domain().communicator().Reduce(&nAtomLocal, &nAtomAll, 1, 
                                        MPI::INT, MPI::SUM, 0);
         #else
         nAtomAll = nAtomLocal;
         #endif
         if (domain().isMaster()) {
            if (nAtomAll != nAtom) {
               UTIL_THROW("nAtomAll != nAtom after distribution");
            }
         }
      }

      int nBond = readGroups<2>(file, "BONDS", "nBond", bondDistributor());

      // Set atom "masks" to suppress pair interactions
      // between covalently bonded atoms.
      if (maskPolicy == MaskBonded) {
         setAtomMasks();
      }

   }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int ConfigIo::writeGroups(std::ostream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& storage,
                  GroupCollector<N>& collector) 
   {
      Group<2>* groupPtr;
      int       nGroup;
      storage.computeNTotal(domain().communicator());
      nGroup = storage.nTotal();
      if (domain().isMaster()) {  
         file << std::endl;
         file << sectionLabel << std::endl;
         file << nGroupLabel << Int(nGroup, 10) << std::endl;
         collector.setup();
         groupPtr = collector.nextPtr();
         while (groupPtr) {
            file << *groupPtr << std::endl;
            groupPtr = collector.nextPtr();
         }
      } else { 
         collector.send();
      }
      return nGroup;
   }

   /* 
   * Write the configuration file.
   */
   void ConfigIo::writeConfig(std::ostream& file)
   {

      // Write Boundary dimensions
      if (domain().isMaster()) {
         file << "BOUNDARY" << std::endl << std::endl;
         file << boundary() << std::endl;
         file << std::endl;
      }

      // Atoms
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  
         file << "ATOMS" << std::endl;
         file << "nAtom" << Int(atomStorage().nAtomTotal(), 10) << std::endl;
         atomCollector_.setup();
         Atom* atomPtr = atomCollector_.nextPtr();
         while (atomPtr) {
            file << Int(atomPtr->id(), 10) << Int(atomPtr->typeId(), 6)
                 << atomPtr->position() << std::endl
                 << "                " << atomPtr->velocity();
            #if 0
            for (int j=0; j < atomPtr->mask().size(); ++j) {
               file << Int(atomPtr->mask()[j], 8);
            }
            #endif
            file << std::endl;
            atomPtr = atomCollector_.nextPtr();
         }
      } else { 
         atomCollector_.send();
      }

      // Write the bonds
      writeGroups<2>(file, "BONDS", "nBond", bondStorage(), bondCollector_);

   }
 
   void ConfigIo::setAtomMasks() 
   {

      AtomIterator     atomIter;
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->mask().clear();
      }
  
      int   atomId0, atomId1; 
      Atom* atomPtr0;
      Atom* atomPtr1;
      GroupIterator<2> bondIter;
      bondStorage().begin(bondIter);
      for ( ; bondIter.notEnd(); ++bondIter) {
         atomId0  = bondIter->atomId(0);
         atomId1  = bondIter->atomId(1);
         atomPtr0 = atomStorage().find(atomId0);
         atomPtr1 = atomStorage().find(atomId1);
         if (atomPtr0) {
            atomPtr0->mask().append(atomId1);
         }
         if (atomPtr1) {
            atomPtr1->mask().append(atomId0);
         }
      }
   }

}
#endif
