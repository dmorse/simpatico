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
      atomStoragePtr_(0)
      #ifdef INTER_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef INTER_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {  
      setClassName("TrajectoryWriter"); 
      domainPtr_ = &simulation.domain();
      boundaryPtr_ = &simulation.boundary();
      Buffer* bufferPtr = &simulation.buffer();

      atomStoragePtr_ = &simulation.atomStorage();
      #ifdef INTER_BOND
      bondStoragePtr_ = &simulation.bondStorage();
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_ = &simulation.angleStorage();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &simulation.dihedralStorage();
      #endif
   }

   /*
   * Read interval and outputFileName. 
   */
   void TrajectoryWriter::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
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

   #if 0
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
   #endif

}
#endif 
