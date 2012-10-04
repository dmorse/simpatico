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

      int nAtom;
      int nBond;
      int nAngle;
      int nDihedral;
      int nImproper;
      int nAtomType;
      int nBondType;
      int nAngleType;
      int nDihedralType;
      int nImproperType;

      if (domain().isMaster()) {  
      
         // Read and discard title line
         std::string       line;
         std::getline(file, line);

         // Read numbers of atoms, bonds, etc.
         file >> nAtom >> Label("atoms");
         file >> nBond >> Label("bonds");
         file >> nAngle >> Label("angles");
         file >> nDihedral >> Label("dihedrals");
         file >> nImproper >> Label("impropers");

         /*
         * Validate nAtom and nBond
         * Lammps files can be read only if the number of atoms and bonds
         * in the lammps file exactly matches the corresponding capacities.
         */
         if (nAtom != atomStorage().totalAtomCapacity()) {
            UTIL_THROW("nAtom != atomCapacity");
         }
         
         if (nBond != bondStorage().totalCapacity()) {
            UTIL_THROW("nBond != bondCapacity");
         }
         

         // Read numbers of atom types, bond types, etc.

         file >> nAtomType >> Label("atom") >> Label("types");
         file >> nBondType >> Label("bond") >> Label("types");
         file >> nAngleType >> Label("angle") >> Label("types");
         file >> nDihedralType >> Label("dihedral") >> Label("types");
         file >> nImproperType >> Label("improper") >> Label("types");

         if (nAtomType > nAtomType_) {
            UTIL_THROW("nAtomType > simulation().nAtomType()");
         }
      }

      // Read and broadcast boundary
      Vector lengths;
      Vector min;
      Vector max;
      if (domain().isMaster()) {  
         file >> min[0] >> max[0] >> Label("xlo") >> Label("xhi");
         file >> min[1] >> max[1] >> Label("ylo") >> Label("yhi");
         file >> min[2] >> max[2] >> Label("zlo") >> Label("zhi");
         lengths.subtract(max, min);
         boundary().setOrthorhombic(lengths);
      }
      #if UTIL_MPI
      bcast(domain().communicator(), boundary(), 0);
      #endif

      // Read particle masses (discard values)
      double mass;
      int    atomTypeId;
      if (domain().isMaster()) {  
         file >> Label("Masses");
         for (int i = 0; i < nAtomType; ++i) {
            file >> atomTypeId >> mass;
         }
      }

      /*
      * Read atomic positions
      *
      * Atom tags must appear in order, numbered from 1
      */
      if (domain().isMaster()) {  
         file >> Label("Atoms");

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
         int moleculeId;
         int rank;
         IntVector shift;

         for(int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            file >> id >> moleculeId >> typeId;
            if (id <= 0 || id > totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id-1);
            atomPtr->setTypeId(typeId-1);
            file >> r;
            atomPtr->position() += min;
            if (UTIL_ORTHOGONAL) {
               atomPtr->position() = r;
            } else {
               boundary().transformCartToGen(r, atomPtr->position());
            }
            file >> shift;

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
         //Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }
        

       
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
   void LammpsOrderedConfigIo::writeConfig(std::ostream& file)
   {
      using std::endl;

      // Atoms
      atomStorage().computeNAtomTotal(domain().communicator());

      // Bonds
      if (bondStorage().capacity()) {
         bondStorage().computeNTotal(domain().communicator());
      }

      if (domain().isMaster()) {
      
         // Write first line (skipped) and a blank line.
         file << "LAMMPS data file" << endl;
         file << endl;
      
         int nAtom = atomStorage().nAtomTotal();

         // Write numbers of atoms, bonds, etc.
         file << Int(nAtom, 10) << " atoms    " << endl;
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
         file << endl;

         IoAtom atom;
         atoms_.reserve(nAtom);
         atoms_.clear();
         atoms_.insert(atoms_.end(), nAtom, atom);

         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         int id;
         int n = 0;
         Vector r;
         int shift;
         while (atomPtr) {
            id = atomPtr->id();
            if (UTIL_ORTHOGONAL) {
               atoms_[id].position = atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
               atoms_[id].position = r;
            }
            atoms_[id].typeId = atomPtr->typeId();
            atoms_[id].id = id;
            atomPtr = atomCollector().nextPtr();
            ++n;
         }
         if (n != nAtom) {
            UTIL_THROW("Something is rotten in Denmark");
         }
         for (id = 0; id < nAtom; ++id) {
            if (id != atoms_[id].id) {
               UTIL_THROW("Something is rotten in Denmark");
            }
            file << Int(id+1, 10) << Int(0,6) << Int(atoms_[id].typeId + 1, 6)
                 << atoms_[id].position;
            for (int i = 0; i < Dimension; ++i) {
               file << Int(shift, 4);
            }
            file << std::endl;
         }
      } else {
         atomCollector().send();
      }

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
