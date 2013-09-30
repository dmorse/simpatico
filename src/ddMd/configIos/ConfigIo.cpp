#ifndef DDMD_CONFIG_IO_CPP
#define DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

#include <ddMd/simulation/Simulation.h>                 
#include <ddMd/communicate/Domain.h>   

#include <ddMd/storage/AtomStorage.h>               
#ifdef INTER_BOND
#include <ddMd/storage/BondStorage.h>               
#endif
#ifdef INTER_ANGLE
#include <ddMd/storage/AngleStorage.h>               
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>               
#endif

#include <ddMd/communicate/Buffer.h> 
#include <ddMd/communicate/GroupCollector.tpp> 
#include <ddMd/communicate/GroupDistributor.tpp> 
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Bond.h>
#include <ddMd/chemistry/MaskPolicy.h>
#include <util/space/Vector.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiLoader.h>
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
      #ifdef INTER_BOND
      bondStoragePtr_(0),
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_(0),
      #endif
      atomCacheCapacity_(0)
      #ifdef INTER_BOND
      , bondCacheCapacity_(0)
      #endif
      #ifdef INTER_ANGLE
      , angleCacheCapacity_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralCacheCapacity_(0)
      #endif
   {  setClassName("ConfigIo"); }

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(Simulation& simulation)
    : domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0),
      #ifdef INTER_BOND
      bondStoragePtr_(0),
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_(0),
      #endif
      atomCacheCapacity_(0)
      #ifdef INTER_BOND
      #endif
      , bondCacheCapacity_(0)
      #ifdef INTER_ANGLE
      , angleCacheCapacity_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralCacheCapacity_(0)
      #endif
   {
      setClassName("ConfigIo"); 
      associate(simulation.domain(),
                simulation.boundary(),
                simulation.atomStorage(),
                #ifdef INTER_BOND
                simulation.bondStorage(),
                #endif
                #ifdef INTER_ANGLE
                simulation.angleStorage(),
                #endif
                #ifdef INTER_DIHEDRAL
                simulation.dihedralStorage(),
                #endif
                simulation.buffer()
               );
   }

   /*
   * Associate with required domain, boundary, storage, and buffer objects.
   */
   void ConfigIo::associate(Domain& domain, Boundary& boundary,
                            AtomStorage& atomStorage,
                            #ifdef INTER_BOND
                            BondStorage& bondStorage,
                            #endif
                            #ifdef INTER_ANGLE
                            AngleStorage& angleStorage,
                            #endif
                            #ifdef INTER_DIHEDRAL
                            DihedralStorage& dihedralStorage,
                            #endif
                            Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;

      atomStoragePtr_ = &atomStorage;
      atomDistributor_.associate(domain, boundary, atomStorage, buffer);
      atomCollector_.associate(domain, atomStorage, buffer);

      #ifdef INTER_BOND
      bondStoragePtr_ = &bondStorage;
      bondDistributor_.associate(domain, atomStorage,
                                 bondStorage, buffer);
      bondCollector_.associate(domain, bondStorage, buffer);
      #endif

      #ifdef INTER_ANGLE
      angleStoragePtr_ = &angleStorage;
      angleDistributor_.associate(domain, atomStorage,
                                  angleStorage, buffer);
      angleCollector_.associate(domain, angleStorage, buffer);
      #endif

      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &dihedralStorage;
      dihedralDistributor_.associate(domain, atomStorage,
                                     dihedralStorage, buffer);
      dihedralCollector_.associate(domain, dihedralStorage, buffer);
      #endif

   }

   /*
   * Set parameters and allocate memory.
   */
   void ConfigIo::initialize(int atomCacheCapacity
                             #ifdef INTER_BOND
                             , int bondCacheCapacity
                             #endif
                             #ifdef INTER_ANGLE
                             , int angleCacheCapacity
                             #endif
                             #ifdef INTER_DIHEDRAL
                             , int dihedralCacheCapacity
                             #endif
                            )
   {
      atomCacheCapacity_ = atomCacheCapacity;
      atomDistributor_.initialize(atomCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);
      #ifdef INTER_BOND
      bondCacheCapacity_ = bondCacheCapacity;
      bondDistributor_.initialize(bondCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
      #endif
      #ifdef INTER_ANGLE
      angleCacheCapacity_ = angleCacheCapacity;
      angleDistributor_.initialize(angleCacheCapacity_);
      angleCollector_.allocate(angleCacheCapacity);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralCacheCapacity_ = dihedralCacheCapacity;
      dihedralDistributor_.initialize(dihedralCacheCapacity_);
      dihedralCollector_.allocate(dihedralCacheCapacity_);
      #endif
   }

   /*
   * Read cache capacity parameters and allocate memory.
   */
   void ConfigIo::readParameters(std::istream& in)
   {
      read<int>(in, "atomCacheCapacity", atomCacheCapacity_);
      atomDistributor_.initialize(atomCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);
      #ifdef INTER_BOND
      read<int>(in, "bondCacheCapacity", bondCacheCapacity_);
      bondDistributor_.initialize(bondCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
      #endif
      #ifdef INTER_ANGLE
      read<int>(in, "angleCacheCapacity", angleCacheCapacity_);
      angleDistributor_.initialize(angleCacheCapacity_);
      angleCollector_.allocate(angleCacheCapacity_);
      #endif
      #ifdef INTER_DIHEDRAL
      read<int>(in, "dihedralCacheCapacity", dihedralCacheCapacity_);
      dihedralDistributor_.initialize(dihedralCacheCapacity_);
      dihedralCollector_.allocate(dihedralCacheCapacity_);
      #endif
   }

   /*
   * Load internal state from input archive and allocate memory.
   */
   void ConfigIo::load(Serializable::IArchive& ar)
   {
      MpiLoader<Serializable::IArchive> loader(*this, ar);

      loader.load(atomCacheCapacity_);
      atomDistributor_.initialize(atomCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);

      #ifdef INTER_BOND
      loader.load(bondCacheCapacity_);
      bondDistributor_.initialize(bondCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
      #endif

      #ifdef INTER_ANGLE
      loader.load(angleCacheCapacity_);
      angleDistributor_.initialize(angleCacheCapacity_);
      angleCollector_.allocate(angleCacheCapacity_);
      #endif

      #ifdef INTER_DIHEDRAL
      loader.load(dihedralCacheCapacity_);
      dihedralDistributor_.initialize(dihedralCacheCapacity_);
      dihedralCollector_.allocate(dihedralCacheCapacity_);
      #endif
   }

   /*
   * Save internal state to output archive.
   */
   void ConfigIo::save(Serializable::OArchive& ar)
   {
      ar & atomCacheCapacity_;
      #ifdef INTER_BOND
      ar & bondCacheCapacity_;
      #endif
      #ifdef INTER_ANGLE
      ar & angleCacheCapacity_;
      #endif
      #ifdef INTER_DIHEDRAL
      ar & dihedralCacheCapacity_;
      #endif
   } 

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int ConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupDistributor<N>& distributor) 
   {
      int nGroup;  // Total number of groups in file
      if (domain().isMaster()) {  
         file >> Label(sectionLabel);
         file >> Label(nGroupLabel) >> nGroup;
         Group<N>* groupPtr;
         int i, j, k;
         distributor.setup();
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
   * Private method to write Group<N> objects.
   */
   template <int N>
   int ConfigIo::writeGroups(std::ofstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& storage,
                  GroupCollector<N>& collector) 
   {
      Group<N>* groupPtr;
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
   * Set Mask (exclusion list) for every atoms, based on bond data.
   */ 
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
