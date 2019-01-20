/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpConfigIo.h"

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
#include <util/misc/FlagSet.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SmpConfigIo::SmpConfigIo(bool hasMolecules)
    : ConfigIo(),
      hasMolecules_(hasMolecules)
   {  setClassName("SmpConfigIo"); }

   /*
   * Constructor.
   */
   SmpConfigIo::SmpConfigIo(Simulation& simulation, bool hasMolecules)
    : ConfigIo(simulation),
      hasMolecules_(hasMolecules)
   {  setClassName("SmpConfigIo"); }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   void SmpConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupDistributor<N>& distributor) 
   {
      if (domain().isMaster()) {  
         int nGroup = 0;  // Total number of groups in file
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
   }

   /*
   * Read a configuration file.
   */
   void SmpConfigIo::readConfig(std::ifstream& file, MaskPolicy maskPolicy)
   {
      // Precondition
      if (atomStorage().nAtom()) {
         UTIL_THROW("Atom storage is not empty (has local atoms)");
      }
      if (atomStorage().nGhost()) {
         UTIL_THROW("Atom storage is not empty (has ghost atoms)");
      }
      if (atomStorage().isCartesian()) {
         UTIL_THROW("Error: Atom storage is set for Cartesian coordinates");
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
      int nAtom = 0;  // Total number of atoms in file
      if (domain().isMaster()) {  

         UTIL_CHECK(Label::isClear());

         // Read ATOMs block header
         file >> Label("ATOMS");

         // Optionally read "ordered" flag
         Label orderedLabel("ordered", false); // optional label
         bool isOrdered = orderedLabel.match(file);

         // Read format string and nAtom
         std::string formatString;
         file >> Label("format") >> formatString;
         UTIL_CHECK(formatString.size() > 0);
   
         // Parse atom format string
         FlagSet atomFormat("imtpvs");
         atomFormat.setActualOrdered(formatString);
         bool hasAtomIndex = atomFormat.isActive('i');
         bool hasAtomContext = atomFormat.isActive('m');
         bool hasAtomTypeId = atomFormat.isActive('t');
         bool hasAtomPosition = atomFormat.isActive('p');
         bool hasAtomVelocity = atomFormat.isActive('v');
         //bool hasAtomShift = atomFormat.isActive('s');
         UTIL_CHECK(hasAtomIndex);
         UTIL_CHECK(hasAtomTypeId);
         UTIL_CHECK(hasAtomPosition);
         UTIL_CHECK(isOrdered);
         // TODO: Add ability to read unordered atoms with context.
         // if (!isOrdered) {
         //   UTIL_CHECK(hasAtomContext);
         // }
         if (hasAtomContext) hasMolecules_ = true;


         // Read number of atoms
         file >> Label("nAtom") >> nAtom;
         // UTIL_CHECK(nAtom == nAtomTot);

         int totalAtomCapacity = atomStorage().totalAtomCapacity();

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().setup();
         #endif

         // Read atoms
         Vector r;
         Vector v;
         Atom*  atomPtr;
         int id;
         int aId;
         int mId;
         int sId;
         int typeId;

         v.zero();
         for (int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();
            
            file >> id;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id);
            if (hasMolecules_) {
               file >> sId >> mId >> aId;
               if (aId < 0) {
                  UTIL_THROW("Invalid Atom");
               }
               if (mId < 0) {
                  UTIL_THROW("Invalid Molecule");
               }
               if (sId < 0) {
                  UTIL_THROW("Invalid Species");
               }
               atomPtr->context().atomId = aId;
               atomPtr->context().moleculeId = mId;
               atomPtr->context().speciesId = sId;
            }
            file >> typeId;
            atomPtr->setTypeId(typeId);
            file >> r;
            boundary().transformCartToGen(r, atomPtr->position());
            if (hasAtomVelocity) {
               file >> v;
            }
            atomPtr->velocity() = v;

            // Add atom to list for sending.
            atomDistributor().addAtom();

         }

         // Send any atoms not sent previously.
         atomDistributor().send();

      } else { // If I am not the master processor
         atomDistributor().receive();
      }

      // Validate atom distribution
      // Checks that all are accounted for and on correct processor
      int nAtomAll;
      nAtomAll = atomDistributor().validate();
      if (domain().isMaster()) {
         if (nAtomAll != nAtom) {
            // nAtom is number of atoms read from file
            // nAtomAll is number found after distribution
            UTIL_THROW("nAtomAll != nAtom after distribution");
         }
      }

      // Read Covalent Groups
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
   int SmpConfigIo::writeGroups(std::ofstream& file, 
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
         collector.setup();
         groupPtr = collector.nextPtr();
         while (groupPtr) {
            file << *groupPtr << std::endl;
            groupPtr = collector.nextPtr();
         }
      } else { 
         collector.send();
      }
      return nGroup;
   }

   /* 
   * Write the configuration file.
   */
   void SmpConfigIo::writeConfig(std::ofstream& file)
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

         file << "ATOMS" << std::endl;

         // Write format string "i[m]tpv"
         std::string format = "format  i";
         if (hasMolecules_) format += "m";
         format += "tpv";
         file << "format" << std::endl;
         file << "nAtom" << Int(atomStorage().nAtomTotal(), 10) << std::endl;

         // Collect and write atoms
         Vector r;
         bool isCartesian = atomStorage().isCartesian();
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            file << Int(atomPtr->id(), 10);
            if (hasMolecules_) {
               file << Int(atomPtr->context().speciesId, 6) 
                    << Int(atomPtr->context().moleculeId, 10)
                    << Int(atomPtr->context().atomId, 6);
            }
            file << Int(atomPtr->typeId(), 6);
            if (isCartesian) {
               r = atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
            }
            file << "\n" << r 
                 << "\n" << atomPtr->velocity() << "\n";
            atomPtr = atomCollector().nextPtr();
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
