#ifndef SPAN_LAMMPS_DUMP_READER_CPP
#define SPAN_LAMMPS_DUMP_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpReader.h" 
#include <spAn/storage/Configuration.h>
#include <util/space/Vector.h>
#include <util/misc/ioUtil.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpReader::LammpsDumpReader(Configuration& configuration)
    : TrajectoryReader(configuration)
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
      line >> step;
      checkString(line, "TIMESTEP");
 
      // Read ITEM: NUMBER OF ATOMS
      notEnd = getNextLine(file, line);
      checkString(line, "ITEM:");
      checkString(line, "NUMBER");
      checkString(line, "OF");
      checkString(line, "ATOMS");
      int nAtom;
      line >> nAtom;

      // Read ITEM: BOX
      notEnd = getNextLine(file, line);
      checkString(line, "ITEM:");
      checkString(line, "BOX");
      // Ignore rest of line
      Vector lengths;
      file >> lengths[0];
      file >> lengths[1];
      file >> lengths[2];
     
      // Read ITEM: ATOMS 
      notEnd = getNextLine(file, line);
      checkString(line, "ITEM:");
      checkString(line, "ATOMS");
      // Ignore the rest of the format line, for now

      // Loop over atoms, read positions
      AtomStorage* storagePtr = &configuration().atoms();
      Atom* atomPtr;
      Vector r;
      IntVector shift;
      int i, j, id, typeId, molId;
      for (i = 0; i < nAtom; ++i) {

         // Find atom address in AtomStorage by id
         file >> id;
         id = id - 1; // Convert from Lammps -> Simpatico convention
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

   }

}
#endif
