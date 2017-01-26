/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpWriter.h"
#include <tools/storage/Configuration.h>
#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Group.h>

#include <util/archives/BinaryFileOArchive.h>
#include <util/format/Dbl.h>
#include <util/space/Vector.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpWriter::LammpsDumpWriter(Processor& processor)
    : TrajectoryWriter(processor, false)
   {}

   /*
   * Constructor.
   */
   LammpsDumpWriter::LammpsDumpWriter(Configuration& configuration,
                                      FileMaster& fileMaster)
    : TrajectoryWriter(configuration, fileMaster, false)
   {}

   /*
   * Destructor.
   */
   LammpsDumpWriter::~LammpsDumpWriter()
   {}

   /*
   *  Write a configuration snapshot. 
   */
   void LammpsDumpWriter::writeFrame(std::ofstream& file, long iStep)
   {

      // Compute total number of atoms
      nAtom_ = atoms().size();

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

      file << "ITEM: ATOMS id type mol x y z" << "\n";
      AtomStorage::Iterator iter;
      atoms().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         id = iter->id;
         typeId = iter->typeId;
         file << id + 1<< " ";
         file << typeId + 1 << " ";
         file << molId << " ";
         for (int i=0; i < Util::Dimension; ++i) {
            file << Dbl(iter->position[i], 13) << " ";
         }
         for (int i=0; i < Util::Dimension; ++i) {
            file << shift << " ";
         }
         file << "\n";
      }

   }

}
