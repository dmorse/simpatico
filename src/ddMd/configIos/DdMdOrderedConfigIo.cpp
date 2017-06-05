/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdOrderedConfigIo.h"

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
   DdMdOrderedConfigIo::DdMdOrderedConfigIo(bool hasMolecules)
    : ConfigIo(),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdOrderedConfigIo"); }

   /*
   * Constructor.
   */
   DdMdOrderedConfigIo::DdMdOrderedConfigIo(Simulation& simulation, bool hasMolecules)
    : ConfigIo(simulation),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdOrderedConfigIo"); }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   void DdMdOrderedConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupDistributor<N>& distributor) 
   {
      if (domain().isMaster()) {  
         int nGroup;  // Total number of groups in file
         file >> Label(sectionLabel);
         file >> Label(nGroupLabel) >> nGroup;

         Group<N>* groupPtr;
         distributor.setup();
         for (int i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            file >> *groupPtr;
            distributor.add();
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         // Receive all groups into BondStorage
         distributor.receive();
      }
      // return nGroup;
   }

   /*
   * Read a configuration file.
   */
   void DdMdOrderedConfigIo::readConfig(std::ifstream& file, MaskPolicy maskPolicy)
   {
      // Precondition
      if (atomStorage().nAtom()) {
         UTIL_THROW("Atom storage is not empty (has local atoms)");
      }
      if (atomStorage().nGhost()) {
         UTIL_THROW("Atom storage is not empty (has ghost atoms)");
      }
      if (atomStorage().isCartesian()) {
         UTIL_THROW("Atom storage is set for Cartesian coordinates");
      }
      if (domain().isMaster() && !file.is_open()) {  
            UTIL_THROW("Error: File is not open on master"); 
      }
      if (!Atom::hasAtomContext()) {
         hasMolecules_ = false;
      }

      // Read and broadcast boundary
      if (domain().isMaster()) {  
         file >> Label("BOUNDARY");
         file >> boundary();
      }
      #if UTIL_MPI
      bcast(domain().communicator(), boundary(), 0);
      #endif

      // Atoms 
      int nAtom;  // Total number of atoms in file
      if (domain().isMaster()) {  

         // Read and distribute atoms
         file >> Label("ATOMS");

         // Read number of atoms
         file >> Label("nAtom") >> nAtom;

         int totalAtomCapacity = atomStorage().totalAtomCapacity();

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().setup();
         #endif

         // Read atoms
         Vector r;
         Atom* atomPtr;
         int id;
         int typeId;

         int aId;
         int mId;
         int sId;

         for (int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();
            
            file >> id >> typeId;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
 
            if (hasMolecules_) {
               file >> sId >> mId >> aId;
               if (sId < 0) {
                  UTIL_THROW("species Id < 0");
               }
               if (mId < 0) {
                  UTIL_THROW("molecule Id < 0");
               }
               if (aId < 0) {
                  UTIL_THROW("atom Id < 0");
               }
               atomPtr->context().speciesId = sId;
               atomPtr->context().moleculeId = mId;
               atomPtr->context().atomId = aId;
            }
  
            file >> r;
            boundary().transformCartToGen(r, atomPtr->position());
            file >> atomPtr->velocity();

            // Add atom to list for sending.
            atomDistributor().addAtom();

         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor
         atomDistributor().receive();
      }

      // Validate atom distribution
      // Check that all atoms are accounted for and on correct processor
      {
         int nAtomAll = atomDistributor().validate();
         if (domain().isMaster()) {
            if (nAtomAll != nAtom) {
               UTIL_THROW("nAtomAll != nAtom after distribution");
            }
         }
      }

      // Read covalent groups
      bool hasGhosts = false;
      #ifdef SIMP_BOND
      if (bondStorage().capacity()) {
         readGroups<2>(file, "BONDS", "nBond", bondDistributor());
         bondStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
         // Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleStorage().capacity()) {
         readGroups<3>(file, "ANGLES", "nAngle", angleDistributor());
         angleStorage().isValid(atomStorage(), domain().communicator(), 
                                hasGhosts);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralStorage().capacity()) {
         readGroups<4>(file, "DIHEDRALS", "nDihedral", dihedralDistributor());
         dihedralStorage().isValid(atomStorage(), domain().communicator(), 
                                   hasGhosts);
      }
      #endif

   }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int DdMdOrderedConfigIo::writeGroups(std::ofstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& storage,
                  GroupCollector<N>& collector) 
   {
      Group<N>* groupPtr;
      int       nGroup;
      storage.computeNTotal(domain().communicator());
      nGroup = storage.nTotal();
      if (domain().isMaster()) {  
         file << std::endl;
         file << sectionLabel << std::endl;
         file << nGroupLabel << Int(nGroup, 10) << std::endl;

         IoGroup<N> ioGroup;
         std::vector<IoGroup <N> > groups;
         groups.reserve(nGroup);
         groups.clear();
         groups.insert(groups.end(), nGroup, ioGroup);

         collector.setup();
         groupPtr = collector.nextPtr();
         int id;
         int n = 0;
         while (groupPtr) {
            id = groupPtr->id();
            groups[id].id = id;
            groups[id].group = *groupPtr;
            groupPtr = collector.nextPtr();
            ++n;
         }
         if (n != nGroup) {
            UTIL_THROW("Something is rotten in Denmark");
         }
         for (id = 0; id < nGroup; ++id) {
            if (id != groups[id].id) {
               UTIL_THROW("Something is rotten in Denmark");
            }
            file << groups[id].group << std::endl;
         }
         file << std::endl;
      } else { 
         collector.send();
      }
      return nGroup;
   }

   /* 
   * Write the configuration file.
   */
   void DdMdOrderedConfigIo::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (domain().isMaster() && !file.is_open()) {  
            UTIL_THROW("Error: File is not open on master"); 
      }
      if (!Atom::hasAtomContext()) {
         hasMolecules_ = false;
      }

      // Write Boundary dimensions
      if (domain().isMaster()) {
         file << "BOUNDARY" << std::endl << std::endl;
         file << boundary() << std::endl;
         file << std::endl;
      }

      // Atoms
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) { 
         int nAtom = atomStorage().nAtomTotal();
         atomCollector().setup();

         file << "ATOMS" << std::endl;
         file << "nAtom" << Int(nAtom, 10) << std::endl;

         // Set up array of nAtom default-constructed elements.
         IoAtom atom;
         atoms_.reserve(nAtom);
         atoms_.clear();
         atoms_.insert(atoms_.end(), nAtom, atom);

         // Collect unordered atoms and store in order in atoms_ .
         Vector r;
         int id;
         int n = 0;
         bool isCartesian = atomStorage().isCartesian();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            if (isCartesian) {
               r = atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
            }
            atoms_[id].position = r;
            atoms_[id].velocity = atomPtr->velocity();
            atoms_[id].typeId = atomPtr->typeId();
            atoms_[id].id = id;
            if (hasMolecules_) {
               atoms_[id].context = atomPtr->context();
            }
            atomPtr = atomCollector().nextPtr();
            ++n;
         } 
         if (n != nAtom) {
            UTIL_THROW("Something is rotten in Denmark");
         }

         // Iterate through atoms_ array to write atoms data.
         for (id = 0; id < nAtom; ++id) {
            if (id != atoms_[id].id) {
               UTIL_THROW("Something is rotten in Denmark");
            }
            file << Int(id, 10) << Int(atoms_[id].typeId, 6);
            if (hasMolecules_) {
               file << Int(atoms_[id].context.speciesId, 6)
                    << Int(atoms_[id].context.moleculeId, 6)
                    << Int(atoms_[id].context.atomId, 6);
            }
            file << "\n" << atoms_[id].position 
                 << "\n" << atoms_[id].velocity 
                 << "\n";
         }

      } else { 
         atomCollector().send();
      }

      // Write the groups
      #ifdef SIMP_BOND
      if (bondStorage().capacity()) {
         writeGroups<2>(file, "BONDS", "nBond", bondStorage(), bondCollector());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleStorage().capacity()) {
         writeGroups<3>(file, "ANGLES", "nAngle", angleStorage(), angleCollector());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralStorage().capacity()) {
         writeGroups<4>(file, "DIHEDRALS", "nDihedral", dihedralStorage(), dihedralCollector());
      }
      #endif

   }
 
}
