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
#include <ddMd/storage/AtomStorage.h>               
#include <ddMd/communicate/Domain.h>   
#include <ddMd/communicate/Buffer.h> 
#include <ddMd/chemistry/Atom.h>
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
      systemPtr_   = &system;
      storagePtr_  = &system.atomStorage();
      domainPtr_   = &system.domain();
      boundaryPtr_ = &system.boundary();
      distributor_.associate(*boundaryPtr_, *domainPtr_, buffer);
   }

   /*
   * Read cache size and allocate memory.
   */
   void ConfigIo::readParam(std::istream& in)
   {
      readBegin(in, "ConfigIo");
      read<int>(in, "cacheCapacity", cacheCapacity_);
      readEnd(in);
      distributor_.setParam(cacheCapacity_);
   }

   /*
   * Read parameters, allocate memory and initialize.
   */
   void ConfigIo::readConfig(std::string filename)
   {

      Atom* ptr;

      int myRank = domain().gridRank();
      int nAtom  = 0;  // Total number of atoms in file

      // If I am the master processor.
      if (myRank == 0) {

         std::ifstream file;

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
         distributor().initSendBuffer();
         #endif

         // Fill the atom objects
         int id;
         int typeId = 0;
         for(int i = 0; i < nAtom; ++i) {

            ptr = distributor().newAtomPtr();

            file >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);
            file >> ptr->position();
            file >> ptr->velocity();

            #if 0
            std::cout << Int(id,6);
            std::cout << Int(typeId,4);
            for (int j = 0; j < Dimension; ++j) {
                std::cout << Dbl(ptr->position()[j], 15, 7);
            }
            for (int j = 0; j < Dimension; ++j) {
                std::cout << Dbl(ptr->velocity()[j], 15, 7);
            }
            std::cout << std::endl;
            #endif

            // Add atom to cache for sending.
            distributor().addAtom(atomStorage());

         }
         file.close();

         // Send any atoms not sent previously.
         distributor().send();

      } else { // If I am not the master processor

         #if UTIL_MPI
         // Receive broadcast of boundary
         bcast(domainPtr_->communicator(), boundary(), 0);
         #endif

         // Receive all atoms into AtomStorage
         distributor().receive(atomStorage());

      }

      // Check that all atoms are accounted for after distribution.
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
