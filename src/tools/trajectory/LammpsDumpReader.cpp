#ifndef TOOLS_LAMMPS_DUMP_READER_CPP
#define TOOLS_LAMMPS_DUMP_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpReader.h" 
#include <tools/storage/Configuration.h>
#include <util/space/Vector.h>
#include <util/misc/ioUtil.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpReader::LammpsDumpReader(Configuration& configuration)
    : TrajectoryReader(configuration, false)
   {  setClassName("LammpsDumpReader"); }

   /*
   * Destructor.
   */
   LammpsDumpReader::~LammpsDumpReader()
   {}

   /*
   * Read a frame.
   */
   bool LammpsDumpReader::readFrame(std::ifstream& file)
   {
      bool notEnd;
      std::stringstream line;

      // Read ITEM: TIMESTEP
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         return false;
      }
      checkString(line, "ITEM:");
      checkString(line, "TIMESTEP");
      int step;
      file >> step;
 
      // Read ITEM: NUMBER OF ATOMS
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: NUMBER OF ATOMS");
      }
      checkString(line, "ITEM:");
      checkString(line, "NUMBER");
      checkString(line, "OF");
      checkString(line, "ATOMS");
      int nAtom;
      file >> nAtom;

      // Read ITEM: BOX
      notEnd = getNextLine(file, line);
      if (!notEnd) {
         UTIL_THROW("EOF reading ITEM: BOX");
      }
      checkString(line, "ITEM:");
      checkString(line, "BOX");
      // Ignore rest of ITEM: BOX line
      Vector min, max, lengths;
      for (int i = 0; i < Dimension; ++i) {
         file >> min[i] >> max[i];
         lengths[i] = max[i] - min[i];
      }

      // Read ITEM: ATOMS 
      notEnd = getNextLine(file, line);
      checkString(line, "ITEM:");
      checkString(line, "ATOMS");
      // Ignore the rest of ITEM: ATOMS  line, for now

      // Loop over atoms, read positions
      AtomStorage* storagePtr = &configuration().atoms();
      Atom* atomPtr;
      Vector r;
      IntVector shift;
      int i, j, id, typeId, molId;
      for (i = 0; i < nAtom; ++i) {

         // Find atom address in AtomStorage by id
         file >> id;
         id = id - 1; // Convert from Lammps -> Simpatico integer convention
         atomPtr = storagePtr->ptr(id);
         if (atomPtr == 0) {
            UTIL_THROW("Unknown atom");
         }
         file >> typeId;
         file >> molId;

         // Read position
         for (j = 0; j < Util::Dimension; ++j) {
            file >> atomPtr->position[j];
         }
         file >> shift;

      }
      return true;
   }

}
#endif
