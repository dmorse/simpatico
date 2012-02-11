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
   ConfigIo::ConfigIo()
    : domainPtr_(0),
      boundaryPtr_(0),
      atomStoragePtr_(0),
      bondStoragePtr_(0),
      atomCacheCapacity_(0),
      bondCacheCapacity_(0)
   {}

   /*
   * Read cache size and allocate memory.
   */
   void ConfigIo::associate(Domain& domain, Boundary& boundary,
                            AtomStorage& atomStorage,
                            BondStorage& bondStorage,
                            Buffer& buffer)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      atomStoragePtr_ = &atomStorage;
      bondStoragePtr_ = &bondStorage;
      atomDistributor_.associate(domain, boundary, buffer);
      bondDistributor_.associate(domain, atomStorage,
                                 bondStorage, buffer);
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

      if (myRank == 0) {
         file.open(filename.c_str());
      }

      // Read and broadcast boundary
      if (myRank == 0) {
         file >> Label("BOUNDARY");
         file >> boundary();
      }
      #if UTIL_MPI
      bcast(domainPtr_->communicator(), boundary(), 0);
      #endif

      // Read and distribute atoms 
      int nAtom;  // Total number of atoms in file
      if (myRank == 0) {

         // Read and distribute atoms
         file >> Label("ATOMS");

         // Read number of atoms
         file >> Label("nAtom") >> nAtom;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " 
                   << nAtom << std::endl;

         int totalAtomCapacity = atomStoragePtr_->totalAtomCapacity();

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().initSendBuffer();
         #endif

         // Read atoms
         Atom* atomPtr;
         int id;
         int typeId;
         int rank;
         for(int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            file >> id >> typeId;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
            file >> atomPtr->position();
            file >> atomPtr->velocity();

            // Add atom to list for sending.
            rank = atomDistributor().addAtom(atomStorage());

         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor

         #if UTIL_MPI
         // Receive all atoms into AtomStorage
         atomDistributor().receive(atomStorage());
         #endif

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
      int nBond;  // Total number of bonds in file
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
         int   i, j, k;
         for (i = 0; i < nBond; ++i) {

            bondPtr = bondDistributor().newPtr();
            file >> *bondPtr;
            for (j = 0; j < 2; ++j) {
               k = bondPtr->atomId(j);
            }
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
         //file << "nAtom" << system().nAtomTotal() << endl;

      }

   }
 
}
#endif
