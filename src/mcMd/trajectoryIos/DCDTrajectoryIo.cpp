#ifndef MCMD_DCD_TRAJECTORY_IO_CPP
#define MCMD_DCD_TRAJECTORY_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DCDTrajectoryIo.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
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
   DCDTrajectoryIo::DCDTrajectoryIo(System &system)
   : TrajectoryIo(system)
   {}

   /*
   * Destructor.
   */
   DCDTrajectoryIo::~DCDTrajectoryIo()
   {}

   static unsigned int read_int(std::fstream &file)
   {
      unsigned int val;
      file.read((char *)&val, sizeof(unsigned int));
      return val;
   }

   void DCDTrajectoryIo::readHeader(std::fstream &file)
   {
      // Calculate atomCapacity for entire simulation
      int atomCapacity = 0;
      int bondCapacity = 0;
      int nSpecies = simulation().nSpecies();
      int speciesCapacity = 0;
      int iSpec,iMol;
      Species* speciesPtr;
      Molecule* molPtr;

      for (iSpec = 0; iSpec < nSpecies; ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         speciesCapacity = speciesPtr->capacity();

         // Add molecules to system
         for (iMol = 0; iMol < speciesCapacity; ++iMol) {
            molPtr = &(speciesPtr->reservoir().pop());
            system().addMolecule(*molPtr);
         }
         atomCapacity += speciesCapacity*speciesPtr->nAtom();
         bondCapacity += speciesCapacity*speciesPtr->nBond();
      }

      // read number of frames
      file.seekp(NFILE_POS);
      nFrames_ = read_int(file);

      file.seekp(NATOMS_POS);
      nAtoms_ = read_int(file);
      if (nAtoms_ != atomCapacity) {
         std::ostringstream oss;
         oss << "Number of atoms in DCD file (" << nAtoms_ << ") does not "
              << "match allocated number of atoms (" << atomCapacity  << ")!";
         UTIL_THROW(oss.str().c_str());
      }

      file.seekp(FRAMEDATA_POS);

      xBuffer_.allocate(nAtoms_);
      yBuffer_.allocate(nAtoms_);
      zBuffer_.allocate(nAtoms_);
   }

   void DCDTrajectoryIo::readFrame(std::fstream &file)
   {
      double lx,ly,lz;
      double angle0,angle1,angle2;

      Vector lengths;

      // Read frame header
      int headerSize;
      headerSize=read_int(file);
      if (headerSize != 6*sizeof(double)) {
         std::ostringstream oss;
         oss << "Unknown file format!";
         UTIL_THROW(oss.str().c_str());
      }

      file.read((char *)&lx, sizeof(double));
      file.read((char *)&angle0, sizeof(double));
      file.read((char *)&ly, sizeof(double));
      file.read((char *)&angle1, sizeof(double));
      file.read((char *)&angle2, sizeof(double));
      file.read((char *)&lz, sizeof(double));

      headerSize=read_int(file);
      if (headerSize != 6*sizeof(double)) {
         std::ostringstream oss;
         oss << "Unkown file format!";
         UTIL_THROW(oss.str().c_str());
      }

      if (!file.good()) {
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
      blockSize = read_int(file);
      file.read((char *) xBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file); // blockSize

      if (blockSize != (int)sizeof(float)*nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize << ", expected " << nAtoms_ << ")";
         UTIL_THROW(oss.str().c_str());
      }

      blockSize = read_int(file);
      file.read((char *) yBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file); // blockSize

      if (blockSize != (int)sizeof(float)*nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize << ", expected " << nAtoms_ << ")";
         UTIL_THROW(oss.str().c_str());
      }

      blockSize = read_int(file);
      file.read((char *) zBuffer_.cArray(), sizeof(float) * nAtoms_);
      read_int(file); // blockSize

      if (blockSize != (int)sizeof(float)* nAtoms_) {
         std::ostringstream oss;
         oss << "Invalid frame size (got " << blockSize << ", expected " << nAtoms_ << ")";
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
   }
}
#endif
