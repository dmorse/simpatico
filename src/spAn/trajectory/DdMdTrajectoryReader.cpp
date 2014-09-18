#ifndef SPAN_DDMD_TRAJECTORY_READER_CPP
#define SPAN_DDMD_TRAJECTORY_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryReader.h" 
#include <spAn/storage/Configuration.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/space/Vector.h>
#include <util/misc/ioUtil.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdTrajectoryReader::DdMdTrajectoryReader(Configuration& configuration)
    : TrajectoryReader(configuration, true)
   {  setClassName("DdMdTrajectoryReader"); }

   /*
   * Destructor.
   */
   DdMdTrajectoryReader::~DdMdTrajectoryReader()
   {}

   void DdMdTrajectoryReader::readHeader(std::ifstream &file)
   {
      BinaryFileIArchive ar(file);
      ar >> nAtom_;
      //std::cout << nAtom_ << std::endl;
   }

   /*
   * Read a frame.
   */
   bool DdMdTrajectoryReader::readFrame(std::ifstream& file)
   {
      BinaryFileIArchive ar(file);

      // Attempt to read iStep
      long iStep; 
      ar >> iStep;
      //std::cout << iStep << std::endl;
      if (file.eof()) {
         return false;
      }

      // Read boundary dimensions
      ar >> configuration().boundary();

      // Loop over atoms, read positions
      AtomStorage* storagePtr = &configuration().atoms();
      Atom* atomPtr;
      int id, i, j;
      for (i = 0; i < nAtom_; ++i) {
         ar >> id;
         //std::cout << id << std::endl;
         atomPtr = storagePtr->ptr(id);
         if (atomPtr == 0) {
            UTIL_THROW("Unknown atom");
         }
         for (j = 0; j < Dimension; ++j) {
            ar >> atomPtr->position[j];
         }
      }

      return true;
   }

}
#endif
