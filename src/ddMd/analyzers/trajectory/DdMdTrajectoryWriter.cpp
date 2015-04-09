/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/space/Vector.h>

#include <climits>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdTrajectoryWriter::DdMdTrajectoryWriter(Simulation& simulation)
    : TrajectoryWriter(simulation, true)
   {  setClassName("DdMdTrajectoryWriter"); }

   /*
   * Destructor.
   */
   DdMdTrajectoryWriter::~DdMdTrajectoryWriter()
   {}

   void DdMdTrajectoryWriter::writeHeader(std::ofstream &file)
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

      if (domain().isMaster()) {  

         ar << iStep;
         ar << boundary();

         Vector r;
         int id, j;
         unsigned int ir;
         bool isCartesian = atomStorage().isCartesian();

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            ar << id;
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
            }
            // ar << atomPtr->velocity();
            atomPtr = atomCollector().nextPtr();
         }

      } else { 
         atomCollector().send();
      }

   }

}
