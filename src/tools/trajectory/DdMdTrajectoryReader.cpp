/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdTrajectoryReader.h" 
#include <tools/storage/Configuration.h>
#include <simp/boundary/Boundary.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/space/Vector.h>
#include <util/misc/ioUtil.h>

#include <climits>

namespace Tools
{

   using namespace Util;
   using namespace Simp;

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
      long iStep = -1;  
      ar >> iStep;

      // Return false if read failed, indicating end of file.
      if (file.eof()) {
         return false;
      }

      // Read boundary dimensions
      Boundary& boundary = configuration().boundary();
      ar >> boundary;

      // Loop over atoms, read atomic positions
      AtomStorage* storagePtr = &configuration().atoms();
      Atom* atomPtr;
      Vector r;
      double h = 1.0/(double(UINT_MAX) + 1.0);
      int id, i, j;
      unsigned int ir;
      for (i = 0; i < nAtom_; ++i) {
         ar >> id;
         atomPtr = storagePtr->ptr(id);
         if (atomPtr == 0) {
            UTIL_THROW("Unknown atom");
         }
         for (j = 0; j < Dimension; ++j) {
            ar >> ir;
            r[j] = ir*h;
         }
         boundary.transformGenToCart(r, atomPtr->position);
      }

      return true;
   }

}
