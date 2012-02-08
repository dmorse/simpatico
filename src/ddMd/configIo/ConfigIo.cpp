#ifndef CONFIG_IO_CPP
#define CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigIo.h"

#include <ddMd/system/System.h>                 
#include <ddMd/communicate/Domain.h>   
#include <ddMd/storage/AtomStorage.h>               
#include <ddMd/storage/BondStorage.h>               
#include <ddMd/communicate/Buffer.h> 
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Bond.h>
#include <util/space/Vector.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigIo::ConfigIo(System& system, Buffer& buffer)
   {
      systemPtr_ = &system;
      domainPtr_ = &system.domain();
      boundaryPtr_ = &system.boundary();
      atomStoragePtr_ = &system.atomStorage();
      bondStoragePtr_ = &system.bondStorage();
      atomDistributor_.associate(*boundaryPtr_, *domainPtr_, buffer);
      bondDistributor_.associate(*bondStoragePtr_, *atomStoragePtr_, 
                                 *domainPtr_, buffer);
      atomCacheCapacity_ = 0;
      bondCacheCapacity_ = 0;
   }

   /*
   * Read cache size and allocate memory.
   */
   void ConfigIo::readParam(std::istream& in)
   {
      readBegin(in, "ConfigIo");
      read<int>(in, "atomCacheCapacity", atomCacheCapacity_);
      read<int>(in, "bondCacheCapacity", bondCacheCapacity_);
      readEnd(in);
      atomDistributor_.setParam(atomCacheCapacity_);
      bondDistributor_.setParam(bondCacheCapacity_);
   }

   /*
   * Read parameters, allocate memory and initialize.
   */
   void ConfigIo::readConfig(std::string filename)
   {

      std::ifstream file;
      int myRank = domain().gridRank();
      int nAtom  = 0;  // Total number of atoms in file
      int nBond  = 0;  // Total number of bonds in file

      // Read and distribute atoms
      if (myRank == 0) {

         file.open(filename.c_str());

         // Read and broadcast system Boundary 
         file >> Label("BOUNDARY");
         file >> boundary();
         #if UTIL_MPI
         bcast(domainPtr_->communicator(), boundary(), 0);
         #endif

         // Read and distribute atoms
         file >> Label("ATOMS");

         // Read number of atoms
         file >> Label("nAtom") >> nAtom;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " 
                   << nAtom << std::endl;

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().initSendBuffer();
         #endif

         // Read atoms
         Atom* atomPtr;
         int id;
         int typeId = 0;
         for(int i = 0; i < nAtom; ++i) {

            atomPtr = atomDistributor().newAtomPtr();

            file >> id >> typeId;
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
            file >> atomPtr->position();
            file >> atomPtr->velocity();

            #if 0
            std::cout << Int(id,6);
            std::cout << Int(typeId,4);
            for (int j = 0; j < Dimension; ++j) {
                std::cout << Dbl(atomPtr->position()[j], 15, 7);
            }
            for (int j = 0; j < Dimension; ++j) {
                std::cout << Dbl(atomPtr->velocity()[j], 15, 7);
            }
            std::cout << std::endl;
            #endif

            // Add atom to cache for sending.
            atomDistributor().addAtom(atomStorage());

         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor

         #if UTIL_MPI
         // Receive broadcast of boundary
         bcast(domainPtr_->communicator(), boundary(), 0);
         #endif

         // Receive all atoms into AtomStorage
         atomDistributor().receive(atomStorage());

      }

      // Check that all atoms are accounted for after distribution.
      {
         int nAtomLocal = atomStorage().nAtom();
         int nAtomAll;
         #ifdef UTIL_MPI
         domain().communicator().Reduce(&nAtomLocal, &nAtomAll, 1, 
                                        MPI::INT, MPI::SUM, 0);
         #else
         nAtomAll = nAtomLocal;
         #endif
         if (myRank == 0) {
            if (nAtomAll != nAtom) {
               UTIL_THROW("nAtomAll != nAtom after distribution");
            }
         }
      }

      // Read and distribute bonds
      if (myRank == 0) {

         // Read and distribute bonds
         file >> Label("BONDS");

         // Read number of bonds
         file >> Label("nBond") >> nBond;

         std::cout << std::endl;
         std::cout << "Num Bonds to be distributed = " 
                   << nBond << std::endl;

         #if UTIL_MPI
         //Initialize the send buffer.
         bondDistributor().initSendBuffer();
         #endif

         // Fill the bond objects
         Bond* bondPtr;
         for (int i = 0; i < nBond; ++i) {

            bondPtr = bondDistributor().newPtr();
            file >> *bondPtr;
            bondDistributor().add();

         }

         // Send any bonds not sent previously.
         bondDistributor().send();

      } else { // If I am not the master processor

         // Receive all bonds into BondStorage
         bondDistributor().receive();

      }

      if (myRank == 0) {
         file.close();
      }

   }

   /* 
   * Write the configuration file.
   */
   void ConfigIo::writeConfig(std::string filename)
   {
      using std::endl;

      int myRank = domain().gridRank();

      // If I am the master processor.
      if (myRank == 0) {

         std::ofstream file;

         // Write Boundary dimensions
         file << "BOUNDARY" << endl << endl;
         file << boundary() << endl;
   
         file << endl << "ATOMS" << endl;
         file << "nAtom" << system().nAtomTotal() << endl;

      }

   }
 
}
#endif
