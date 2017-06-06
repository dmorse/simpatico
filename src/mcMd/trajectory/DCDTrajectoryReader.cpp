/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DCDTrajectoryReader.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <vector>
#include <sstream>

//! File position of NFILE in DCD header
#define NFILE_POS 8L
//! File position of NATOMS in DCD header
#define NATOMS_POS 268L
//! File position for start of frame data
#define FRAMEDATA_POS 276L

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   DCDTrajectoryReader::DCDTrajectoryReader(System &system)
   : TrajectoryReader(system),
     nAtoms_(0),
     nFrames_(0),
     frameId_(0)
   {}

   /*
   * Destructor.
   */
   DCDTrajectoryReader::~DCDTrajectoryReader()
   {}

   static unsigned int read_int(std::fstream &file)
   {
      unsigned int val;
      file.read((char *)&val, sizeof(unsigned int));
      return val;
   }

   void DCDTrajectoryReader::open(std::string filename)
   {
      // Open trajectory file
      file_.open(filename.c_str(), std::ios::in | std::ios::binary);
      if (file_.fail()) {
         std::string message;
         message= "Error opening trajectory file. Filename: " + filename;
         UTIL_THROW(message.c_str());
      }

      // Read number of frames
      file_.seekp(NFILE_POS);
      nFrames_ = read_int(file_);
      frameId_ = 0;

      // Read number of atoms
      file_.seekp(NATOMS_POS);
      nAtoms_ = read_int(file_);

      // Add all molecules and check consistency
      addMolecules();

      if (nAtoms_ != nAtomTotal_) {
         std::ostringstream oss;
         oss << "Number of atoms in DCD file (" << nAtoms_ << ") does not "
              << "match allocated number of atoms (" << nAtomTotal_  << ")!";
         UTIL_THROW(oss.str().c_str());
      }

      file_.seekp(FRAMEDATA_POS);

      xBuffer_.allocate(nAtoms_);
      yBuffer_.allocate(nAtoms_);
      zBuffer_.allocate(nAtoms_);
   }

   bool DCDTrajectoryReader::readFrame()
   {
      // Check if the last frame was already read
      if (frameId_ >= nFrames_) {
         return false;
      }

      double lx, ly, lz;
      double angle0, angle1, angle2;

      Vector lengths;

      // Read frame header
      int headerSize;
      headerSize=read_int(file_);
      if (headerSize != 6*sizeof(double)) {
         std::ostringstream oss;
         oss << "Unknown file format!";
         UTIL_THROW(oss.str().c_str());
      }

      file_.read((char *)&lx, sizeof(double));
      file_.read((char *)&angle0, sizeof(double));
      file_.read((char *)&ly, sizeof(double));
      file_.read((char *)&angle1, sizeof(double));
      file_.read((char *)&angle2, sizeof(double));
      file_.read((char *)&lz, sizeof(double));

      headerSize=read_int(file_);
      if (headerSize != 6*sizeof(double)) {
         std::ostringstream oss;
         oss << "Unkown file format!";
         UTIL_THROW(oss.str().c_str());
      }

      if (!file_.good()) {
         std::ostringstream oss;
         oss << "Error reading trajectory file!";
         UTIL_THROW(oss.str().c_str());
      }

      lengths[0]=lx;
      lengths[1]=ly;
      lengths[2]=lz;
      boundary().setOrthorhombic(lengths);

      // Read frame data
      int blockSize;

      // read coords
      blockSize = read_int(file_);
      file_.read((char *) xBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file_); // blockSize

      if (blockSize != (int)sizeof(float)*nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize 
             << ", expected " << nAtoms_ << ")";
         UTIL_THROW(oss.str().c_str());
      }

      blockSize = read_int(file_);
      file_.read((char *) yBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file_); // blockSize

      if (blockSize != (int)sizeof(float)*nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize 
             << ", expected " << nAtoms_ << ")";
         UTIL_THROW(oss.str().c_str());
      }

      blockSize = read_int(file_);
      file_.read((char *) zBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file_); // blockSize

      if (blockSize != (int)sizeof(float)* nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize 
             << ", expected " << nAtoms_ << ")";
         UTIL_THROW(oss.str().c_str());
      }

      // Load positions, assume they are ordered according to species
      int iSpecies,iMol;
      int bufferIdx=0;
      Species *speciesPtr;
      Molecule::AtomIterator atomIter;
      Molecule *molPtr;

      for (iSpecies = 0; iSpecies < simulation().nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         for (iMol = 0; iMol < speciesPtr->capacity(); ++iMol) {
            molPtr = &system().molecule(iSpecies, iMol);
            for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->position()[0] = (double) xBuffer_[bufferIdx];
               atomIter->position()[1] = (double) yBuffer_[bufferIdx];
               atomIter->position()[2] = (double) zBuffer_[bufferIdx];

               // shift into simulation cell
               boundary().shift(atomIter->position());

               bufferIdx++;
            }
         }
      }

      // Increment frame counter
      ++frameId_;

      // Indicate successful completion
      return true;
   }

   void DCDTrajectoryReader::close()
   { file_.close(); }

}
