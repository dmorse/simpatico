/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SerializeConfigIo.h"

#include <ddMd/simulation/Simulation.h>                 
#include <ddMd/communicate/Domain.h>   

#include <ddMd/storage/AtomStorage.h>               
#ifdef SIMP_BOND
#include <ddMd/storage/BondStorage.h>               
#endif
#ifdef SIMP_ANGLE
#include <ddMd/storage/AngleStorage.h>               
#endif
#ifdef SIMP_DIHEDRAL
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
   * Private method to load Group<N> objects. (Call on all processors).
   */
   template <int N>
   int SerializeConfigIo::loadGroups(Serializable::IArchive& ar,
                                     GroupDistributor<N>& distributor) 
   {
      int nGroup = 0;  // Total number of groups in archive
      if (domain().isMaster()) {  
         ar >> nGroup;
         Group<N>* groupPtr;
         distributor.setup();
         for (int i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            ar >> *groupPtr;
            distributor.add();
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         // Receive all groups into BondStorage
         distributor.receive();
      }
      return nGroup; // Valid only on master
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
      if (atomStorage().isCartesian()) {
         UTIL_THROW("Atom storage set for Cartesian coordinates");
      }

      // Load and broadcast boundary
      if (domain().isMaster()) {  
         ar >> boundary();
      }
      #if UTIL_MPI
      bcast(domain().communicator(), boundary(), 0);
      #endif

      // Load atoms 
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
         Atom* atomPtr;
         AtomContext* contextPtr;
         int  id;
         int  typeId;
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
            ar >> atomPtr->groups();
            if (Atom::hasAtomContext()) {
               contextPtr = &atomPtr->context();
               ar >> contextPtr->speciesId;
               ar >> contextPtr->moleculeId;
               ar >> contextPtr->atomId;
            }
            ar >> r;
            boundary().transformCartToGen(r, atomPtr->position());
            ar >> atomPtr->velocity();

            // Add atom to list for sending.
            atomDistributor().addAtom();
         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor
         atomDistributor().receive();
      }

      // Validate atom distribution:
      // Check that all are accounted for and on correct processor
      int nAtomAll;
      nAtomAll = atomDistributor().validate();
      if (domain().isMaster()) {
         if (nAtomAll != nAtom) {
            UTIL_THROW("nAtomAll != nAtom after distribution");
         }
      }

      // Load groups
      bool hasGhosts = false;
      #ifdef SIMP_BOND
      if (bondStorage().capacity()) {
         loadGroups<2>(ar, bondDistributor());
         bondStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
         // Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleStorage().capacity()) {
         loadGroups<3>(ar, angleDistributor());
         angleStorage().isValid(atomStorage(), domain().communicator(), 
                                hasGhosts);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
   * Save the configuration to an archive.
   */
   void SerializeConfigIo::saveConfig(Serializable::OArchive& ar)
   {
      // Write Boundary dimensions
      if (domain().isMaster()) {
         ar << boundary();
      }

      // Save atoms
      bool isCartesian = atomStorage().isCartesian();
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {  

         int id;
         int typeId;
         int nAtom = atomStorage().nAtomTotal();
         Vector r;
         AtomContext* contextPtr;

         ar << nAtom;
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            typeId = atomPtr->typeId();
            ar << id;
            ar << typeId;
            ar << atomPtr->groups();
            if (Atom::hasAtomContext()) {
               contextPtr = &atomPtr->context();
               ar << contextPtr->speciesId;
               ar << contextPtr->moleculeId;
               ar << contextPtr->atomId;
            }
            if (isCartesian) {
               ar << atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
               ar << r;
            }
            ar << atomPtr->velocity();
            atomPtr = atomCollector().nextPtr();
         }

      } else { 
         atomCollector().send();
      }

      // Save groups
      #ifdef SIMP_BOND
      if (bondStorage().capacity()) {
         saveGroups<2>(ar, bondStorage(), bondCollector());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleStorage().capacity()) {
         saveGroups<3>(ar, angleStorage(), angleCollector());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralStorage().capacity()) {
         saveGroups<4>(ar, dihedralStorage(), dihedralCollector());
      }
      #endif
   }
 
   /*
   * Read configuration file, using an input archive.
   */
   void SerializeConfigIo::readConfig(std::ifstream& file, MaskPolicy maskPolicy)
   {
      // Precondition
      if (domain().isMaster() && !file.is_open()) {  
            UTIL_THROW("Error: File is not open on master"); 
      }
      // Other preconditions are enforced by loadConfig
      
      Serializable::IArchive ar(file);
      loadConfig(ar, maskPolicy);
   }

   /*
   * Write configuration file, using an output file archive.
   */
   void SerializeConfigIo::writeConfig(std::ofstream& file)
   {
      // Preconditions
      if (domain().isMaster() && !file.is_open()) {  
            UTIL_THROW("Error: File is not open on master"); 
      }
      
      Serializable::OArchive ar(file);
      saveConfig(ar);
   }

}
