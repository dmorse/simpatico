/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdGroupTrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/space/Vector.h>
#include <util/misc/Bit.h>

#include <climits>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdGroupTrajectoryWriter::DdMdGroupTrajectoryWriter(Simulation& simulation)
    : TrajectoryWriter(simulation, true)
   {}

   /*
   * Destructor.
   */
   DdMdGroupTrajectoryWriter::~DdMdGroupTrajectoryWriter()
   {}

   /*
   * Read parameter file block.
   */
   void DdMdGroupTrajectoryWriter::readParameters(std::istream& in)
   {
      TrajectoryWriter::readParameters(in);
      read<unsigned int>(in, "groupId", groupId_);
   }

   /*
   * Load internal state from an archive.
   */
   void DdMdGroupTrajectoryWriter::loadParameters(Serializable::IArchive &ar)
   {
      TrajectoryWriter::loadParameters(ar);
      loadParameter<unsigned int>(ar, "groupId", groupId_);
   }

   /*
   * Save internal state to output archive.
   */
   void DdMdGroupTrajectoryWriter::save(Serializable::OArchive& ar)
   {
      TrajectoryWriter::save(ar);
      ar << groupId_;
   }

   void DdMdGroupTrajectoryWriter::writeHeader(std::ofstream &file)
   {
      if (domain().isMaster()) {  
         BinaryFileOArchive ar(file);
         Bit bit(groupId_);

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         nAtom_ = 0;
         while (atomPtr) {
            if (bit.isSet(atomPtr->groups())) {
               ++nAtom_;
            }
            atomPtr = atomCollector().nextPtr();
         }
         ar << nAtom_;
         //file << nAtom_ << std::endl;
      } else { 
         atomCollector().send();
      }
   }

   void DdMdGroupTrajectoryWriter::writeFrame(std::ofstream &file, long iStep)
   {
      if (domain().isMaster()) {  
         BinaryFileOArchive ar(file);
         Bit bit(groupId_);

         ar << iStep;
         ar << boundary();
         //file << iStep << std::endl;
         //file << boundary() << std::endl;

         Vector r;
         int id, j;
         unsigned int ir;
         bool isCartesian = atomStorage().isCartesian();

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            if (bit.isSet(atomPtr->groups())) {
               id = atomPtr->id();
               ar << id;
               //file << id << std::endl;
               if (isCartesian) {
                  boundary().transformCartToGen(atomPtr->position(), r);
               } else {
                  r = atomPtr->position();
               }
               for (j = 0; j < Dimension; ++j) {
                  if (r[j] >= 1.0) r[j] -= 1.0;
                  if (r[j] <  0.0) r[j] += 1.0;
                  ir = floor( UINT_MAX*r[j] + r[j] + 0.5 );
                  ar << ir;
                  //file << ir;
               }
               // ar << atomPtr->velocity();
            }
            atomPtr = atomCollector().nextPtr();
         }
      } else {
         atomCollector().send();
      }
   }

}
