#ifndef DDMD_DDMD_TRAJECTORY_WRITER_CPP
#define DDMD_DDMD_TRAJECTORY_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <vector>
#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdTrajectoryWriter::DdMdTrajectoryWriter(Simulation& simulation)
    : TrajectoryWriter(simulation)
   {}

   /*
   * Destructor.
   */
   DdMdTrajectoryWriter::~DdMdTrajectoryWriter()
   {}

   void DdMdTrajectoryWriter::writeHeader(std::ofstream &file, long iStep)
   {
      BinaryFileOArchive ar(file);

      // Compute and write total number of atoms
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  
         nAtom_ = atomStorage().nAtomTotal();
         ar << nAtom_;
      }

      // Write groups
      #ifdef INTER_BOND
      if (bondStorage().capacity()) {
         writeGroups<2>(ar, bondStorage(), bondCollector());
      }
      #endif
      #ifdef INTER_ANGLE
      if (angleStorage().capacity()) {
         writeGroups<3>(ar, angleStorage(), angleCollector());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralStorage().capacity()) {
         writeGroups<4>(ar, dihedralStorage(), dihedralCollector());
      }
      #endif
   }

   void DdMdTrajectoryWriter::writeFrame(std::ofstream &file, long iStep)
   {
      BinaryFileOArchive ar(file);

      // Write Boundary dimensions
      if (domain().isMaster()) {
         ar << boundary();
      }

      // Writes atoms
      // atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  
         int id;
         int typeId;
         Vector r;
         bool isCartesian = atomStorage().isCartesian();

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            typeId = atomPtr->typeId();
            ar << id;
            ar << typeId;
            if (isCartesian) {
               ar << atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
               ar << r;
            }
            ar << atomPtr->velocity();
            atomPtr = atomCollector().nextPtr();
         }

      } else { 
         atomCollector().send();
      }

   }

   /*
   * Private method to save Group<N> objects.
   */
   template <int N>
   int DdMdTrajectoryWriter::writeGroups(BinaryFileOArchive& ar,
                  GroupStorage<N>& storage, GroupCollector<N>& collector) 
   {
      Group<N>* groupPtr;
      int       nGroup;
      storage.computeNTotal(domain().communicator());
      nGroup = storage.nTotal();
      if (domain().isMaster()) {  
         ar << nGroup;
         collector.setup();
         groupPtr = collector.nextPtr();
         while (groupPtr) {
            ar << *groupPtr;
            groupPtr = collector.nextPtr();
         }
      } else { 
         collector.send();
      }
      return nGroup;
   }

}
#endif
