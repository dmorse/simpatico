/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

#include <ddMd/simulation/Simulation.h>                 
#include <ddMd/communicate/Domain.h>   
#include <ddMd/communicate/GroupDistributor.tpp>   
#include <ddMd/communicate/GroupCollector.tpp>   

#include <ddMd/storage/AtomStorage.h>               
#ifdef SIMP_BOND
#include <ddMd/storage/BondStorage.h>               
//#include <ddMd/chemistry/Bond.h>
#endif
#ifdef SIMP_ANGLE
#include <ddMd/storage/AngleStorage.h>               
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>               
#endif

#include <ddMd/communicate/Buffer.h> 
#include <ddMd/communicate/GroupCollector.tpp> 
#include <ddMd/communicate/GroupDistributor.tpp> 
#include <ddMd/chemistry/Atom.h>
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
      atomStoragePtr_(0)
      #ifdef SIMP_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {  setClassName("ConfigIo"); }

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(Simulation& simulation)
    : domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0)
      #ifdef SIMP_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {
      setClassName("ConfigIo"); 
      associate(simulation.domain(),
                simulation.boundary(),
                simulation.atomStorage(),
                #ifdef SIMP_BOND
                simulation.bondStorage(),
                #endif
                #ifdef SIMP_ANGLE
                simulation.angleStorage(),
                #endif
                #ifdef SIMP_DIHEDRAL
                simulation.dihedralStorage(),
                #endif
                simulation.buffer()
               );
   }

   /*
   * Destructor.
   */
   ConfigIo::~ConfigIo()
   {}

   /*
   * Associate with required domain, boundary, storage, and buffer objects.
   */
   void ConfigIo::associate(Domain& domain, Boundary& boundary,
                            AtomStorage& atomStorage,
                            #ifdef SIMP_BOND
                            BondStorage& bondStorage,
                            #endif
                            #ifdef SIMP_ANGLE
                            AngleStorage& angleStorage,
                            #endif
                            #ifdef SIMP_DIHEDRAL
                            DihedralStorage& dihedralStorage,
                            #endif
                            Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;

      atomStoragePtr_ = &atomStorage;
      #ifdef SIMP_BOND
      bondStoragePtr_ = &bondStorage;
      #endif
      #ifdef SIMP_ANGLE
      angleStoragePtr_ = &angleStorage;
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralStoragePtr_ = &dihedralStorage;
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
   * Set Mask (exclusion list) for all atoms, based on bond data.
   */ 
   void ConfigIo::setAtomMasks() 
   {

      AtomIterator     atomIter;
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->mask().clear();
      }
  
      #ifdef SIMP_BOND
      int   atomId0, atomId1;
      Atom* atomPtr0;
      Atom* atomPtr1;
      GroupIterator<2> bondIter;
      const AtomMap& atomMap = atomStorage().map();
      bondStorage().begin(bondIter);
      for ( ; bondIter.notEnd(); ++bondIter) {
         atomId0  = bondIter->atomId(0);
         atomId1  = bondIter->atomId(1);
         atomPtr0 = atomMap.find(atomId0);
         atomPtr1 = atomMap.find(atomId1);
         if (atomPtr0) {
            atomPtr0->mask().append(atomId1);
         }
         if (atomPtr1) {
            atomPtr1->mask().append(atomId0);
         }
      }
      #endif
   }

}
