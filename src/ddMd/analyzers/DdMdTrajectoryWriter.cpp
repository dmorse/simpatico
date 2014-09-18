#ifndef DDMD_DDMD_TRAJECTORY_WRITER_CPP
#define DDMD_DDMD_TRAJECTORY_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/space/Vector.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdTrajectoryWriter::DdMdTrajectoryWriter(Simulation& simulation)
    : TrajectoryWriter(simulation, true)
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

   }

   void DdMdTrajectoryWriter::writeFrame(std::ofstream &file, long iStep)
   {
      BinaryFileOArchive ar(file);

      // Write iStep and boundary dimensions
      if (domain().isMaster()) {
         ar << iStep;
         ar << boundary();
      }

      // Write atoms
      if (domain().isMaster()) {  
         Vector r;
         bool isCartesian = atomStorage().isCartesian();

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            ar << atomPtr->id();
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

}
#endif
