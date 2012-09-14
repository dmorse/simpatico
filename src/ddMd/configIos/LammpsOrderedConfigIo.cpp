#ifndef DDMD_LAMMPS_CONFIG_IO_CPP
#define DDMD_LAMMPS_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsOrderedConfigIo.h"

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
   LammpsOrderedConfigIo::LammpsOrderedConfigIo()
    : ConfigIo()
   {}

   /*
   * Constructor.
   */
   LammpsOrderedConfigIo::LammpsOrderedConfigIo(Simulation& simulation)
    : ConfigIo(simulation)
   {
      nAtomType_ =  simulation.nAtomType();
      nBondType_ =  simulation.nBondType();
   }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int LammpsOrderedConfigIo::readGroups(std::istream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupDistributor<N>& distributor) 
   {
      int nGroup;  // Total number of groups in file
      if (domain().isMaster()) {  
         file >> Label(sectionLabel);
         file >> Label(nGroupLabel) >> nGroup;
         Group<N>* groupPtr;
         int i, j, k;
         distributor.setup();
         for (i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            file >> *groupPtr;
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
   * Read a configuration file.
   */
   void LammpsOrderedConfigIo::readConfig(std::istream& file, MaskPolicy maskPolicy)
   {
      // Precondition
      if (atomStorage().nAtom()) {
         UTIL_THROW("Atom storage is not empty (has local atoms)");
      }
      if (atomStorage().nGhost()) {
         UTIL_THROW("Atom storage is not empty (has ghost atoms)");
      }
      if (UTIL_ORTHOGONAL) {
         if (!atomStorage().isCartesian()) {
            UTIL_THROW("Atom storage must use Cartesian coordinates");
         }
      } else {
         if (atomStorage().isCartesian()) {
            UTIL_THROW("Atom storage must use generalized coordinates");
         }
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
         Atom*  atomPtr;
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
            file >> r;
            if (UTIL_ORTHOGONAL) {
               atomPtr->position() = r;
            } else {
               boundary().transformCartToGen(r, atomPtr->position());
            }
            file >> atomPtr->velocity();

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

      bool hasGhosts = false;

      if (bondStorage().capacity()) {
         readGroups<2>(file, "BONDS", "nBond", bondDistributor());
         bondStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
         // Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }

      #ifdef INTER_ANGLE
      if (angleStorage().capacity()) {
         readGroups<3>(file, "ANGLES", "nAngle", angleDistributor());
         angleStorage().isValid(atomStorage(), domain().communicator(), 
                                hasGhosts);
      }
      #endif

      #ifdef INTER_DIHEDRAL
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
   int LammpsOrderedConfigIo::writeGroups(std::ostream& file, 
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
   void LammpsOrderedConfigIo::writeConfig(std::ostream& file)
   {
      using std::endl;

      // Write first line (skipped) and a blank line.
      file << "LAMMPS data file" << endl;
      file << endl;

      // Atoms
      atomStorage().computeNAtomTotal(domain().communicator());

      // Bonds
      if (bondStorage().capacity()) {
         bondStorage().computeNTotal(domain().communicator());
      }

      if (domain().isMaster()) {

      // Write numbers of atoms, bonds, etc.
      file << Int(atomStorage().nAtomTotal(), 10) << " atoms    " << endl;
      file << Int(bondStorage().nTotal(), 10) << " bonds    " << endl;
      file << Int(0)     << " angles   " << endl;
      file << Int(0)     << " dihedrals" << endl;
      file << Int(0)     << " impropers" << endl;
      file << endl;

      // Write numbers of atom types, bond types, etc.
      file << Int(nAtomType_) << " atom types" << endl;
      file << Int(nBondType_) << " bond types" << endl;
      file << Int(0) << " angle types" << endl;
      file << Int(0) << " dihedral types" << endl;
      file << Int(0) << " improper types" << endl;
      file << endl;

      // Write Boundary dimensions
      Vector lengths = boundary().lengths();
      file << Dbl(0.0) << Dbl(lengths[0]) << "  xlo xhi" << endl;
      file << Dbl(0.0) << Dbl(lengths[1]) << "  ylo yhi" << endl;
      file << Dbl(0.0) << Dbl(lengths[2]) << "  zlo zhi" << endl;
      file << endl;

      // Write masses (all set to 1.0 for now)
      // lammps atom type = Simpatico atom type + 1
      file << "Masses" << endl;
      file << endl;
      for (int iType = 0; iType < nAtomType_; ++iType) {
          file << Int(iType+1, 5) << Dbl(1.0) << endl;
      }
      file << endl;

      // Write atomic positions
      // lammps atom     tag = Simpatico atom id + 1
      // lammps molecule id  = Simpatico molecule id + 1
      file << "Atoms" << endl;
       atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         int       shift;
         while (atomPtr) {
            file << Int(atomPtr->id(), 10) << Int(atomPtr->typeId(), 6);
             file     << atomPtr->position();
               for (int i = 0; i < Dimension; ++i) {
                  file << Int(shift, 4);
               }
            file << std::endl;
            atomPtr = atomCollector().nextPtr();
         }
     } else {
         atomCollector().send();
      }


      // Write bond topology
      file << "Bonds" << endl;
      file << endl;
      // Write the groups
      if (bondStorage().capacity()) {
         writeGroups<2>(file, "BONDS", "nBond", bondStorage(), bondCollector());
      }

      /*
      // Write bond topology
      file << "Bonds" << endl;
      file << endl;
      Molecule::BondIterator bondIter;
      int                    iBond = 1;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
               file << Int(iBond, 8 ) << Int(bondIter->typeId() + 1, 5);
               file << Int(bondIter->atom(0).id() + 1, 8);
               file << Int(bondIter->atom(1).id() + 1, 8);
               file << endl;
               ++iBond;
            }
         }
      }
      file << endl;
      */
   }
 
}
#endif
