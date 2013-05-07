#ifndef DDMD_SERIALIZE_CONFIG_IO_CPP
#define DDMD_SERIALIZE_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SerializeConfigIo.h"

#include <ddMd/simulation/Simulation.h>                 
#include <ddMd/communicate/Domain.h>   

#include <ddMd/storage/AtomStorage.h>               
#include <ddMd/storage/BondStorage.h>               
#ifdef INTER_ANGLE
#include <ddMd/storage/AngleStorage.h>               
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>               
#endif

#include <ddMd/communicate/GroupCollector.tpp> 
#include <ddMd/communicate/GroupDistributor.tpp> 

#include <ddMd/communicate/Buffer.h> 
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Bond.h>
#include <ddMd/chemistry/MaskPolicy.h>
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
   SerializeConfigIo::SerializeConfigIo()
    : ConfigIo()
   {  setClassName("SerializeConfigIo"); }

   /*
   * Constructor.
   */
   SerializeConfigIo::SerializeConfigIo(Simulation& simulation)
    : ConfigIo(simulation)
   {  setClassName("SerializeConfigIo"); }

   /*
   * Private method to load Group<N> objects.
   */
   template <int N>
   int SerializeConfigIo::loadGroups(Serializable::IArchive& ar,
                                     GroupDistributor<N>& distributor) 
   {
      int nGroup;  // Total number of groups in archive
      if (domain().isMaster()) {  
         ar >> nGroup;
         Group<N>* groupPtr;
         int i, j, k;
         distributor.setup();
         for (i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            ar >> *groupPtr;
            for (j = 0; j < 2; ++j) {
               k = groupPtr->atomId(j);
            }
            distributor.add();
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         // Receive all groups into BondStorage
         distributor.receive();
      }
      return nGroup;
   }

   /*
   * Load a configuration from input archive.
   */
   void SerializeConfigIo::loadConfig(Serializable::IArchive& ar, MaskPolicy maskPolicy)
   {
      // Preconditions
      if (atomStorage().nAtom()) {
         UTIL_THROW("Atom storage is not empty (has local atoms)");
      }
      if (atomStorage().nGhost()) {
         UTIL_THROW("Atom storage is not empty (has ghost atoms)");
      }
      if (UTIL_ORTHOGONAL) {
         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Atom coordinates are not Cartesian");
         }
      } else {
         if (atomStorage().isCartesian()) {
            UTIL_THROW("Atom coordinates are Cartesian");
         }
      }

      // Read and broadcast boundary
      if (domain().isMaster()) {  
         ar >> boundary();
      }
      #if UTIL_MPI
      bcast(domain().communicator(), boundary(), 0);
      #endif

      // Atoms 
      int nAtom;  // Total number of atoms in archive
      if (domain().isMaster()) {  

         ar >> nAtom;
         int totalAtomCapacity = atomStorage().totalAtomCapacity();

         #ifdef UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().setup();
         #endif

         // Read atoms
         Vector  r;
         Atom*  atomPtr;
         int  id;
         int  typeId;
         int  rank;
         for (int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            ar >> id;
            ar >> typeId;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            if (typeId < 0) {
               UTIL_THROW("Negative atom type id");
            }
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
            ar >> r;
            if (UTIL_ORTHOGONAL) {
               atomPtr->position() = r;
            } else {
               boundary().transformCartToGen(r, atomPtr->position());
            }
            ar >> atomPtr->velocity();

            // Add atom to list for sending.
            rank = atomDistributor().addAtom();
         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor
         atomDistributor().receive();
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
         if (domain().isMaster()) {
            if (nAtomAll != nAtom) {
               UTIL_THROW("nAtomAll != nAtom after distribution");
            }
         }
      }

      // Load groups
      bool hasGhosts = false;
      if (bondStorage().capacity()) {
         loadGroups<2>(ar, bondDistributor());
         bondStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
         // Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }
      #ifdef INTER_ANGLE
      if (angleStorage().capacity()) {
         loadGroups<3>(ar, angleDistributor());
         angleStorage().isValid(atomStorage(), domain().communicator(), 
                                hasGhosts);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralStorage().capacity()) {
         loadGroups<4>(ar, dihedralDistributor());
         dihedralStorage().isValid(atomStorage(), domain().communicator(), 
                                   hasGhosts);
      }
      #endif

   }

   /*
   * Private method to save Group<N> objects.
   */
   template <int N>
   int SerializeConfigIo::saveGroups(Serializable::OArchive& ar,
                  GroupStorage<N>& storage, GroupCollector<N>& collector) 
   {
      Group<N>* groupPtr;
      int       nGroup;
      storage.computeNTotal(domain().communicator());
      nGroup = storage.nTotal();
      if (domain().isMaster()) {  
         ar << nGroup;
         collector.setup();
         groupPtr = collector.nextPtr();
         while (groupPtr) {
            ar << *groupPtr;
            groupPtr = collector.nextPtr();
         }
      } else { 
         collector.send();
      }
      return nGroup;
   }

   /* 
   * Write the configuration to archive.
   */
   void SerializeConfigIo::saveConfig(Serializable::OArchive& ar)
   {
      // Write Boundary dimensions
      if (domain().isMaster()) {
         ar << boundary();
      }

      // Atoms
      bool isCartesian = atomStorage().isCartesian();
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  

         int id;
         int typeId;
         int nAtom = atomStorage().nAtomTotal();
         Vector r;

         ar << nAtom;
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            typeId = atomPtr->typeId();
            ar << id;
            ar << typeId;
            if (isCartesian) {
               ar << atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
               ar << r;
            }
            ar << atomPtr->velocity();
            #if 0
            for (int j=0; j < atomPtr->mask().size(); ++j) {
               ar << atomPtr->mask()
            }
            #endif
            atomPtr = atomCollector().nextPtr();
         }
      } else { 
         atomCollector().send();
      }

      // Write the groups
      if (bondStorage().capacity()) {
         saveGroups<2>(ar, bondStorage(), bondCollector());
      }
      #ifdef INTER_ANGLE
      if (angleStorage().capacity()) {
         saveGroups<3>(ar, angleStorage(), angleCollector());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralStorage().capacity()) {
         saveGroups<4>(ar, dihedralStorage(), dihedralCollector());
      }
      #endif

   }
 
   /*
   * Read configuration file.
   *
   * This routine opens and reads a file on the master, and distributes
   * atom data among the processors.
   */
   void SerializeConfigIo::readConfig(std::ifstream& file, MaskPolicy maskPolicy)
   {
       #if 0
       std::ifstream* ptr = dynamic_cast<std::ifstream*>(&file);
       if (!ptr) {
          UTIL_THROW("Failed dynamic cast: Input istream in not a std::ifstream");
       }
       if (!ptr->is_open()) {
          UTIL_THROW("File not open for reading");
       }
       #endif
       Serializable::IArchive ar(file);
       loadConfig(ar, maskPolicy);
   }

   /*
   * Write configuration file.
   *
   * This routine opens and writes a file on the master,
   * collecting atom data from all processors.
   */
   void SerializeConfigIo::writeConfig(std::ofstream& file)
   {
       #if 0
       std::ofstream* ptr = dynamic_cast<std::ofstream*>(&file);
       if (!ptr) {
          UTIL_THROW("Failed dynamic cast: Output ostream in not a std::ofstream");
       }
       if (!ptr->is_open()) {
          UTIL_THROW("File not open for writing");
       }
       #endif
       Serializable::OArchive ar(file);
       saveConfig(ar);
   }
   
}
#endif
