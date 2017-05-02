/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsConfigIo.h"

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
   LammpsConfigIo::LammpsConfigIo()
    : ConfigIo()
   {}

   /*
   * Constructor.
   */
   LammpsConfigIo::LammpsConfigIo(Simulation& simulation)
    : ConfigIo(simulation)
   {
      nAtomType_ = simulation.nAtomType();
      nBondType_ = 0;
      nAngleType_ = 0;
      nDihedralType_ = 0;
      nImproperType_ = 0;
      #ifdef SIMP_BOND
      nBondType_ = simulation.nBondType();
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_ = simulation.nAngleType();
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_ = simulation.nDihedralType();
      #endif
   }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   void LammpsConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel, int nGroup,
                  GroupDistributor<N>& distributor) 
   {
      if (domain().isMaster()) {  
         file >> Label(sectionLabel);
         Group<N>* groupPtr;
         int i, j, k;
         distributor.setup();
         for (i = 0; i < nGroup; ++i) {
            groupPtr = distributor.newPtr();
            file >> k;
            groupPtr->setId(k-1);
            file >> k;
            groupPtr->setTypeId(k-1);
            for (j = 0; j < N; ++j) {
               file >> k;
               groupPtr->setAtomId(j, k-1);
            }
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
   void LammpsConfigIo::readConfig(std::ifstream& file, MaskPolicy maskPolicy)
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

         // Read numbers of atom types, bond types, etc.
         file >> nAtomType >> Label("atom") >> Label("types");
         file >> nBondType >> Label("bond") >> Label("types");
         file >> nAngleType >> Label("angle") >> Label("types");
         file >> nDihedralType >> Label("dihedral") >> Label("types");
         file >> nImproperType >> Label("improper") >> Label("types");

         if (nAtomType > nAtomType_) {
            UTIL_THROW("nAtomType > simulation().nAtomType()");
         }
         #ifdef SIMP_BOND
         if (nBondType > nBondType_) {
            UTIL_THROW("nAtomType > simulation().nBondType()");
         }
         #endif
         #ifdef SIMP_ANGLE
         if (nAngleType > nAngleType_) {
            UTIL_THROW("nAngleype > simulation().nAngleType()");
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (nDihedralType > nDihedralType_) {
            UTIL_THROW("nDihedralType > simulation().nDihedralType()");
         }
         #endif
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

         //Initialize the send buffer.
         atomDistributor().setup();
         
         // Read atoms
         Vector r(0.0);
         Atom*  atomPtr = 0;
         int id;
         int typeId;
         int moleculeId;
         IntVector shift;

         for (int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            file >> id >> moleculeId >> typeId;
            if (id <= 0 || id > totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            atomPtr->setId(id-1);
            atomPtr->setTypeId(typeId-1);
            file >> r;
            atomPtr->position() += min; // Shift corner of Boundary to (0, 0, 0)
            boundary().transformCartToGen(r, atomPtr->position());
            file >> shift;

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
      int nAtomAll;
      nAtomAll = atomDistributor().validate();
      if (domain().isMaster()) {
         if (nAtomAll != nAtom) {
            UTIL_THROW("nAtomAll != nAtom after distribution");
         }
      }

      bool hasGhosts = false;
      #ifdef SIMP_BOND
      if (bondStorage().capacity()) {
         readGroups<2>(file, "Bonds", nBond, bondDistributor());
         bondStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
         //Set atom "mask" values
         if (maskPolicy == MaskBonded) {
            setAtomMasks();
         }
      }
      #endif
       
      #ifdef SIMP_ANGLE
      if (angleStorage().capacity()) {
         readGroups<3>(file, "Angles", nAngle, angleDistributor());
         angleStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
      }
      #endif
       
      #ifdef SIMP_DIHEDRAL
      if (dihedralStorage().capacity()) {
         readGroups<4>(file, "Dihedrals", nDihedral, dihedralDistributor());
         dihedralStorage().isValid(atomStorage(), domain().communicator(), hasGhosts);
      }
      #endif
       
   }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   void LammpsConfigIo::writeGroups(std::ofstream& file, 
                  const char* sectionLabel,
                  GroupStorage<N>& storage,
                  GroupCollector<N>& collector) 
   {
      Group<N>* groupPtr;
      int       nGroup;
      storage.computeNTotal(domain().communicator());
      nGroup = storage.nTotal();
      if (domain().isMaster()) { 

         IoGroup<N> ioGroup;
         std::vector<IoGroup <N> > groups;

         // Collect and sort groups
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
            UTIL_THROW("Inconsistency in total number of groups");
         }

         // Write groups
         file << std::endl;
         file << sectionLabel << std::endl;
         file << std::endl;
         int j, k;
         for (id = 0; id < nGroup; ++id) {
            if (id != groups[id].id) {
               UTIL_THROW("Incorrect group id in ordered output");
            }
            k = groups[id].id + 1;
            file << k;
            k = groups[id].group.typeId() + 1;
            file << " " << k;
            for (j = 0; j < N; ++j) {
               k = groups[id].group.atomId(j) + 1;
               file << " " << k ;
            }
            file << std::endl;
         }
         file << std::endl;
      } else { 
         collector.send();
      }
   }

   /* 
   * Write the configuration file.
   */
   void LammpsConfigIo::writeConfig(std::ofstream& file)
   {
      // Preconditions
      if (domain().isMaster() && !file.is_open()) {  
            UTIL_THROW("Error: File is not open on master"); 
      }

      using std::endl;

      // Atoms
      atomStorage().computeNAtomTotal(domain().communicator());

      // Bonds
      #ifdef SIMP_BOND
      if (nBondType_) {
         if (bondStorage().capacity()) {
            bondStorage().computeNTotal(domain().communicator());
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         if (angleStorage().capacity()) {
            angleStorage().computeNTotal(domain().communicator());
         }
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         if (dihedralStorage().capacity()) {
            dihedralStorage().computeNTotal(domain().communicator());
         }
      }
      #endif

      if (domain().isMaster()) {
      
         // Write first line (skipped) and a blank line.
         file << "LAMMPS data file" << endl;
         file << endl;
      
         int nAtom = atomStorage().nAtomTotal();
         int nBond = 0;
         int nAngle  = 0;
         int nDihedral= 0;
         int nImproper = 0;
         #ifdef SIMP_BOND
         if (nBondType_) {
            nBond = bondStorage().nTotal();
         }
         #endif
         #ifdef SIMP_ANGLE
         if (nAngleType_) {
            nAngle = angleStorage().nTotal();
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (nDihedralType_) {
            nDihedral = dihedralStorage().nTotal();
         }
         #endif

         // Write numbers of atoms, bonds, etc.
         file << nAtom << " atoms     "      << endl;
         file << nBond << " bonds     " << endl;
         file << nAngle << " angles    " << endl;
         file << nDihedral << " dihedrals " << endl;
         file << nImproper << " impropers" << endl;

         // Write numbers of atom types, bond types, etc.
         file << endl;
         file << nAtomType_     << " atom types" << endl;
         file << nBondType_     << " bond types" << endl;
         file << nAngleType_    << " angle types" << endl;
         file << nDihedralType_ << " dihedral types" << endl;
         file << nImproperType_ << " improper types" << endl;

         // Write Boundary dimensions
         file << endl;
         Vector lengths = boundary().lengths();
         file << Dbl(0.0) << Dbl(lengths[0]) << "  xlo xhi" << endl;
         file << Dbl(0.0) << Dbl(lengths[1]) << "  ylo yhi" << endl;
         file << Dbl(0.0) << Dbl(lengths[2]) << "  zlo zhi" << endl;

         // Write masses (all set to 1.0 for now)
         // lammps atom type = Simpatico atom type + 1
         file << endl;
         file << "Masses" << endl;
         file << endl;
         for (int iType = 0; iType < nAtomType_; ++iType) {
            file << iType+1 << " " << 1.0 << endl;
         }

         IoAtom atom;
         atoms_.reserve(nAtom);
         atoms_.clear();
         atoms_.insert(atoms_.end(), nAtom, atom);


         // Collect atoms
         atomCollector().setup();
         Vector r;
         int id;
         int n = 0;
         bool isCartesian = atomStorage().isCartesian();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            id = atomPtr->id();
            atoms_[id].id = id;
            atoms_[id].typeId = atomPtr->typeId();
            if (isCartesian) {
               r = atomPtr->position();
            } else {
               boundary().transformGenToCart(atomPtr->position(), r);
            }
            atoms_[id].position = r;
            atomPtr = atomCollector().nextPtr();
            ++n;
         }
         if (n != nAtom) {
            UTIL_THROW("Inconsistency in number of atoms");
         }

         // Write atom positions
         // lammps atom     tag = Simpatico atom id + 1
         // lammps molecule id  = Simpatico molecule id + 1
         file << endl;
         file << "Atoms" << endl;
         file << endl;
         int shift = 0;
         for (id = 0; id < nAtom; ++id) {
            if (id != atoms_[id].id) {
               UTIL_THROW("Incorrect atom id in ordered output");
            }
            file << id+1 << " " << "1 " << atoms_[id].typeId + 1 
                 << " " << atoms_[id].position;
            for (int i = 0; i < Dimension; ++i) {
               file << " " << shift;
            }
            file << std::endl;
         }

         // Write atomic velocities
         // lammps atom     tag = Simpatico atom id + 1
         // lammps molecule id  = Simpatico molecule id + 1
         /*
         file << endl;
         file << "Velocities" << endl;
         file << endl;
         for (id = 0; id < nAtom; ++id) {
            if (id != atoms_[id].id) {
               UTIL_THROW("Incorrect atom id in ordered output");
            }
            file << id+1 << " " << atoms_[id].velocity;
            file << std::endl;
         }
         */
      } else {
         atomCollector().send();
      }

      // Write the groups
      #ifdef SIMP_BOND
      if (nBondType_) {
         if (bondStorage().capacity()) {
            writeGroups<2>(file, "Bonds", bondStorage(), bondCollector());
         }
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         if (angleStorage().capacity()) {
            writeGroups<3>(file, "Angles", angleStorage(), angleCollector());
         }
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         if (dihedralStorage().capacity()) {
            writeGroups<4>(file, "Dihedrals", dihedralStorage(), dihedralCollector());
         }
      }
      #endif

   }
 
}
