#ifndef DDMD_TRAJECTORY_WRITER_CPP
#define DDMD_TRAJECTORY_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
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
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Bond.h>

#include <util/space/Vector.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiLoader.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryWriter::TrajectoryWriter(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false),
      domainPtr_(0),
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
   {  
      setClassName("TrajectoryWriter"); 
      domainPtr_ = &simulation.domain();
      boundaryPtr_ = &simulation.boundary();

      atomStoragePtr_ = &simulation.atomStorage();
      atomCollector_.associate(*domainPtr_, *atomStoragePtr_, buffer);

      #ifdef INTER_BOND
      bondStoragePtr_ = &simulation.bondStorage();
      bondCollector_.associate(*domainPtr_, *bondStoragePtr_, buffer);
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_ = &simulation.angleStorage();
      angleCollector_.associate(*domainPtr_, *angleStoragePtr_, buffer);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &simulation.dihedralStorage();
      dihedralCollector_.associate(*domainPtr_, *dihedralStoragePtr_, buffer);
      #endif

   }

   /*
   * Read interval and outputFileName. 
   */
   void TrajectoryWriter::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      bool isRequired = false;
      atomCacheCapacity_ = 500;
      read<int>(in, "atomCacheCapacity", atomCacheCapacity_, isRequired);
      atomCollector_.allocate(atomCacheCapacity_);
      #ifdef INTER_BOND
      bondCacheCapacity_ = 500;
      read<int>(in, "bondCacheCapacity", bondCacheCapacity_, isRequired);
      bondCollector_.allocate(bondCacheCapacity_);
      #endif
      #ifdef INTER_ANGLE
      angleCacheCapacity_ = 500;
      read<int>(in, "angleCacheCapacity", angleCacheCapacity_, isRequired);
      angleCollector_.allocate(angleCacheCapacity_);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralCacheCapacity_ = 500;
      read<int>(in, "dihedralCacheCapacity", dihedralCacheCapacity_, isRequired);
      dihedralCollector_.allocate(dihedralCacheCapacity_);
      #endif
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void TrajectoryWriter::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      loader.load(atomCacheCapacity_);
      atomCollector_.allocate(atomCacheCapacity_);

      #ifdef INTER_BOND
      loader.load(bondCacheCapacity_);
      bondCollector_.allocate(bondCacheCapacity_);
      #endif
      #ifdef INTER_ANGLE
      loader.load(angleCacheCapacity_);
      angleCollector_.allocate(angleCacheCapacity_);
      #endif
      #ifdef INTER_DIHEDRAL
      loader.load(dihedralCacheCapacity_);
      dihedralCollector_.allocate(dihedralCacheCapacity_);
      #endif

      isInitialized_ = true;
   }

   /*
   * Save internal state to output archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive& ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
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
   * Clear sample counter and close file.
   */
   void TrajectoryWriter::clear() 
   {  
      nSample_ = 0; 
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
   }

   /*
   * Write frame to file, header on first sample.
   */
   void TrajectoryWriter::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         if (nSample_ == 0) {
            simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
            writeHeader(outputFile_, iStep);
         }
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }

   /*
   * Clear sample counter and close output file.
   */
   void TrajectoryWriter::output()
   {  clear(); }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int TrajectoryWriter::writeGroups(std::ofstream& file, 
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

}
#endif 
