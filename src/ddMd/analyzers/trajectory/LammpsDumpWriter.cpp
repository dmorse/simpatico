/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/format/Dbl.h>
#include <util/space/Vector.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpWriter::LammpsDumpWriter(Simulation& simulation)
    : TrajectoryWriter(simulation)
   {  setClassName("LammpsDumpWriter"); }

   /*
   * Destructor.
   */
   LammpsDumpWriter::~LammpsDumpWriter()
   {}

   /*
   *  Write a configuration snapshot. 
   */
   void LammpsDumpWriter::writeFrame(std::ofstream &file, long iStep)
   {

      // Compute total number of atoms
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  
         nAtom_ = atomStorage().nAtomTotal();
      }

      if (domain().isMaster()) {
         file << "ITEM: TIMESTEP" << "\n";
         file << iStep << "\n";

         file << "ITEM: NUMBER OF ATOMS" << "\n";
         file << nAtom_ << "\n";

         file << "ITEM: BOX BOUNDS pp pp pp" << "\n";
         Vector lengths = boundary().lengths();
         file << Dbl(0.0) << Dbl(lengths[0]) << "\n";
         file << Dbl(0.0) << Dbl(lengths[1]) << "\n";
         file << Dbl(0.0) << Dbl(lengths[2]) << "\n";

         Vector r;
         int id;
         int typeId;
         int shift = 0;
         int molId = 1;
         bool isCartesian = atomStorage().isCartesian();

         file << "ITEM: ATOMS id type mol x y z" << "\n";
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            typeId = atomPtr->typeId();
            if (isCartesian) {
               r = atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
            }
            file << id + 1<< " ";
            file << typeId + 1 << " ";
            file << molId << " ";
            for (int i=0; i < Util::Dimension; ++i) {
               file << Dbl(atomPtr->position()[i], 13) << " ";
            }
            for (int i=0; i < Util::Dimension; ++i) {
               file << shift << " ";
            }
            file << "\n";
            atomPtr = atomCollector().nextPtr();
         }

      } else { 
         atomCollector().send();
      }

   }

}
