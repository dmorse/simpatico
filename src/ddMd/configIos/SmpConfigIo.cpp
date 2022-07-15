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

#include <simp/species/Species.h>

#include <util/space/Vector.h>
#include <util/param/Label.h>
#include <util/param/OptionalLabel.h>
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
   SmpConfigIo::SmpConfigIo(Simulation& simulation)
    : ConfigIo(simulation),
      simulationPtr_(&simulation)
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
   * Private method to replicate and send Group<N> objects for a species.
   */
   template <int N>
   void
   SmpConfigIo::sendSpeciesGroups(int& groupId, int& firstAtomId,
                                  int  nMolecule, int  nAtom, int  nGroup,
                                  const DArray<SpeciesGroup<N> >& groups,
                                  GroupDistributor<N>& distributor)
   {
      if (domain().isMaster()) {
         Group<N>* groupPtr = 0;
         const SpeciesGroup<N>* speciesGroupPtr = 0;
         int im, ig, j, ia;
         for (im = 0; im < nMolecule; ++im) {
            for (ig = 0; ig < nGroup; ++ig) {
               groupPtr = distributor.newPtr();
               groupPtr->setId(groupId);
               speciesGroupPtr = &groups[ig];
               groupPtr->setTypeId(speciesGroupPtr->typeId());
               for (j = 0; j < N; ++j) {
                  ia = firstAtomId + speciesGroupPtr->atomId(j);
                  groupPtr->setAtomId(j, ia);
               }
               //Log::file() << *groupPtr << std::endl;
               distributor.add();
               ++groupId;
            } // end loop over groups
            firstAtomId += nAtom;
         } // end loop over molecules
      }
   }

   #ifdef SIMP_BOND
   void SmpConfigIo::makeBonds()
   {
      GroupDistributor<2>& distributor = bondDistributor();
      if (domain().isMaster()) {
         int nSpecies = simulation().nSpecies();
         int groupId = 0;
         int firstAtomId = 0;
         distributor.setup();
         for (int is = 0; is < nSpecies; ++is) {
            Species& species = simulation().species(is);
            sendSpeciesGroups(groupId,
                              firstAtomId,
                              species.capacity(),
                              species.nAtom(),
                              species.nBond(),
                              species.speciesBonds(),
                              distributor);
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         distributor.receive();
      }
   }
   #endif

   #ifdef SIMP_ANGLE
   void SmpConfigIo::makeAngles()
   {
      GroupDistributor<3>& distributor = angleDistributor();
      if (domain().isMaster()) {
         int nSpecies = simulation().nSpecies();
         int groupId = 0;
         int firstAtomId = 0;
         distributor.setup();
         for (int is = 0; is < nSpecies; ++is) {
            Species& species = simulation().species(is);
            sendSpeciesGroups(groupId,
                              firstAtomId,
                              species.capacity(),
                              species.nAtom(),
                              species.nAngle(),
                              species.speciesAngles(),
                              distributor);
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         distributor.receive();
      }
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   void SmpConfigIo::makeDihedrals()
   {
      GroupDistributor<4>& distributor = dihedralDistributor();
      if (domain().isMaster()) {
         int nSpecies = simulation().nSpecies();
         int groupId = 0;
         int firstAtomId = 0;
         distributor.setup();
         for (int is = 0; is < nSpecies; ++is) {
            Species& species = simulation().species(is);
            sendSpeciesGroups(groupId,
                              firstAtomId,
                              species.capacity(),
                              species.nAtom(),
                              species.nDihedral(),
                              species.speciesDihedrals(),
                              distributor);
         }
         // Send any groups not sent previously.
         distributor.send();
      } else { // If I am not the master processor
         distributor.receive();
      }
   }
   #endif

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
      if (domain().isMaster()) {
         if (!file.is_open()) {
            UTIL_THROW("Error: File is not open for reading on master");
         }
         if (!Label::isClear()) {
            UTIL_THROW("Error: Label buffer is not clear on master");
         }
      }

      int nSpecies = 0;  // Number of species
      int nAtomTot = 0;  // Number of molecules of all species
      DArray<int> firstAtomIds;

      // Optionally read SPECIES block
      if (domain().isMaster()) {
         OptionalLabel speciesBlockLabel("SPECIES");
         if (speciesBlockLabel.match(file)) {
            file >> Label("nSpecies") >> nSpecies;
            UTIL_CHECK(nSpecies > 0);
         }
         firstAtomIds.allocate(nSpecies);
      }
      #if UTIL_MPI
      bcast(domain().communicator(), nSpecies, 0);
      #endif
      if (nSpecies > 0) {

         // Set nSpecies in parent Simulation.
         // If a Species array was allocated previously, reallocate.
         simulation().setNSpecies(nSpecies);

         // Read molecular structures for all species on master
         if (domain().isMaster()) {
            Label speciesLabel("species");
            Label nMoleculeLabel("nMolecule");
            int j, nMolecule;
            firstAtomIds[0] = 0;
            for (int i = 0; i < nSpecies; ++i) {
               file >> speciesLabel >> j;
               UTIL_CHECK(j == i);
               firstAtomIds[i] = nAtomTot;
               file >> nMoleculeLabel >> nMolecule;
               UTIL_CHECK(nMolecule > 0);
               Species& species = simulation().species(i);
               species.readStructure(file);
               species.setCapacity(nMolecule);
               nAtomTot += nMolecule * species.nAtom();
            }
            // Note: Species information is only stored on the master
            // processor and is not broadcast or replicated. The value
            // of nSpecies, however, is broadcast and stored on all
            // processors.
         }
      }

      // BOUNDARY block
      if (domain().isMaster()) {
         file >> Label("BOUNDARY");
         file >> boundary();
      }
      #if UTIL_MPI
      bcast(domain().communicator(), boundary(), 0);
      #endif

      // ATOMS block
      int nAtom = 0;  // Total number of atoms in file
      if (domain().isMaster()) {

         UTIL_CHECK(Label::isClear());

         // Read ATOMs block header
         file >> Label("ATOMS");

         // Optionally read "ordered" flag
         OptionalLabel orderedLabel("ordered");
         bool isOrdered = orderedLabel.match(file);

         // Default format settings
         bool hasAtomContext = Atom::hasAtomContext();
         bool hasAtomVelocity = true;

         // Optionally read atom format string
         OptionalLabel formatLabel("format");
         bool hasFormat = formatLabel.match(file);
         if (hasFormat) {
            std::string formatString;
            file >> formatString;
            UTIL_CHECK(formatString.size() > 0);

            // Parse and validate atom format string
            FlagSet atomFormat("itmpvs");
            atomFormat.setActualOrdered(formatString);
            bool hasAtomIndex = atomFormat.isActive('i');
            UTIL_CHECK(hasAtomIndex);
            bool hasAtomTypeId = atomFormat.isActive('t');
            UTIL_CHECK(hasAtomTypeId);
            bool hasAtomPosition = atomFormat.isActive('p');
            UTIL_CHECK(hasAtomPosition);
            hasAtomContext = atomFormat.isActive('m');
            UTIL_CHECK(hasAtomContext == Atom::hasAtomContext());
            hasAtomVelocity = atomFormat.isActive('v');
            //bool hasAtomShift = atomFormat.isActive('s');
         }

         // TODO: Add ability to read unordered atoms with context.
         // if (!isOrdered) {
         //   UTIL_CHECK(hasAtomContext);
         // }

         // Read number of atoms
         file >> Label("nAtom") >> nAtom;
         // if (nSpecies) {
         //    UTIL_CHECK(nAtom == nAtomTot);
         // }

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor().setup();
         #endif

         Vector r;
         Vector v;
         Atom*  atomPtr;
         int id;
         int aId;
         int mId;
         int sId;
         int typeId;
         int totalAtomCapacity = atomStorage().totalAtomCapacity();

         // Loop over atoms
         v.zero();
         for (int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor().newAtomPtr();

            // Atom index
            file >> id;
            if (id < 0 || id >= totalAtomCapacity) {
               UTIL_THROW("Invalid atom id");
            }
            if (isOrdered) {
               UTIL_CHECK(id == i);
            }
            atomPtr->setId(id);

            // Atom type id
            file >> typeId;
            atomPtr->setTypeId(typeId);

            // Atom context
            if (hasAtomContext) {
               file >> sId >> mId >> aId;
               if (sId < 0) {
                  UTIL_THROW("Invalid Species");
               }
               if (mId < 0) {
                  UTIL_THROW("Invalid Molecule");
               }
               if (aId < 0) {
                  UTIL_THROW("Invalid Atom");
               }
               atomPtr->context().speciesId = sId;
               atomPtr->context().moleculeId = mId;
               atomPtr->context().atomId = aId;
            }
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

      // Read or replicate covalent groups:
      // If nSpecies > 0, then replicate and broadcast groups.
      // If nSpecies == 0, then read and broadcast all groups.
      bool hasGhosts = false;

      if (simulation().nSpecies() == 0) {

         // Read Covalent Groups
         #ifdef SIMP_BOND
         if (bondStorage().capacity()) {
            readGroups<2>(file, "BONDS", "nBond", bondDistributor());
            bondStorage().isValid(atomStorage(), domain().communicator(),
                                  hasGhosts);
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
            readGroups<4>(file, "DIHEDRALS", "nDihedral",
                          dihedralDistributor());
            dihedralStorage().isValid(atomStorage(), domain().communicator(),
                                      hasGhosts);
         }
         #endif

      } else {     // if simulation().hasSpecies()

         #ifdef SIMP_BOND
         if (bondStorage().capacity()) {
            makeBonds();
            bondStorage().isValid(atomStorage(), domain().communicator(),
                                  hasGhosts);
            // Set atom "mask" values
            if (maskPolicy == MaskBonded) {
               setAtomMasks();
            }
         }
         #endif
         #ifdef SIMP_ANGLE
         if (angleStorage().capacity()) {
            makeAngles();
            angleStorage().isValid(atomStorage(), domain().communicator(),
                                   hasGhosts);
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (dihedralStorage().capacity()) {
            makeDihedrals();
            dihedralStorage().isValid(atomStorage(), domain().communicator(),
                                   hasGhosts);
         }
         #endif

      }

      // Postcondition: Label buffer should be clear on exit
      if (domain().isMaster()) {
         UTIL_CHECK(Label::isClear());
      }
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

      // Compute total total number of groups on all processors
      // nGroup = total only on master, zero on other processors
      storage.computeNTotal(domain().communicator());
      int nGroup = storage.nTotal();

      // Collect groups and write to file
      if (domain().isMaster()) {
         using std::endl;
         file << endl;
         file << sectionLabel << endl;
         file << nGroupLabel << Int(nGroup, 10) << endl;
         Group<N>* groupPtr;
         collector.setup();
         groupPtr = collector.nextPtr();
         while (groupPtr) {
            file << *groupPtr << endl;
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
      using std::endl;

      // Write SPECIES block (if any)
      if (domain().isMaster()) {
         if (simulation().hasSpecies()) {
            file << "SPECIES" << endl << endl;
            int nSpecies = simulation().nSpecies();
            UTIL_CHECK(nSpecies > 0);
            for (int i = 0; i < nSpecies; ++i) {
               file << "species  " << i << endl;
               Species& species = simulation().species(i);
               file << "mNolecule  " << species.capacity() << endl;
               species.writeStructure(file);
               file << endl << endl;
            }
         }
      }

      // Write BOUNDARY block (periodic unit cell dimensions)
      if (domain().isMaster()) {
         file << "BOUNDARY" << endl;
         file << boundary() << endl;
         file << endl;
      }

      // WRITE ATOMS block
      bool hasAtomContext = Atom::hasAtomContext();
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {

         file << "ATOMS" << endl;

         // Note: No "ordered" flag is written, because this function
         // writes unordered atoms.

         // Write format string "it[m]pv"
         std::string format = "format  it";
         if (hasAtomContext) format += "m";
         format += "pv";
         file << format << endl;
         file << "nAtom" << Int(atomStorage().nAtomTotal(), 10)
              << endl;

         // Collect and write atoms
         Vector r;
         bool isCartesian = atomStorage().isCartesian();
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while (atomPtr) {
            file << Int(atomPtr->id(), 10);
            file << Int(atomPtr->typeId(), 6);
            if (hasAtomContext) {
               file << Int(atomPtr->context().speciesId, 6)
                    << Int(atomPtr->context().moleculeId, 10)
                    << Int(atomPtr->context().atomId, 6);
            }
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

      if (!simulation().hasSpecies()) {
         // Write the covalent groups
         #ifdef SIMP_BOND
         if (bondStorage().capacity()) {
            writeGroups<2>(file, "BONDS", "nBond", bondStorage(),
                           bondCollector());
         }
         #endif
         #ifdef SIMP_ANGLE
         if (angleStorage().capacity()) {
            writeGroups<3>(file, "ANGLES", "nAngle", angleStorage(),
                           angleCollector());
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (dihedralStorage().capacity()) {
            writeGroups<4>(file, "DIHEDRALS", "nDihedral", dihedralStorage(),
                           dihedralCollector());
         }
         #endif
      }

   }

}
