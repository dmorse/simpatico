/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulation.h"

// namespace McMd
#include <mcMd/analyzers/AnalyzerManager.h>
#include <mcMd/species/SpeciesManager.h>
#include <mcMd/chemistry/Activate.h>

// namespace Simp
#include <simp/species/Species.h>
#include <simp/species/SpeciesGroup.tpp>

// namespace Util
#include <util/containers/ArraySet.h>
#include <util/param/Factory.h>
#include <util/space/Vector.h>
#include <util/misc/initStatic.h>

#ifdef UTIL_MPI
#include "McMd_mpi.h"
#endif

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   #ifdef UTIL_MPI
   Simulation::Simulation(MPI::Intracomm& communicator)
    : iStep_(0),
      nSystem_(1),
      speciesManagerPtr_(0),
      analyzerManagerPtr_(0),
      moleculeCapacity_(0),
      nAtomType_(-1),
      atomCapacity_(0)
      #ifndef SIMP_NOPAIR
      , maskedPairPolicy_(MaskBonded)
      #endif
      #ifdef SIMP_BOND
      , nBondType_(-1)
      , bondCapacity_(0)
      #endif
      #ifdef SIMP_ANGLE
      , nAngleType_(-1)
      , angleCapacity_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , nDihedralType_(-1)
      , dihedralCapacity_(0)
      #endif
      #ifdef SIMP_COULOMB
      , hasCoulomb_(-1)
      #endif
      #ifdef SIMP_EXTERNAL
      , hasExternal_(-1)
      #endif
      #ifdef MCMD_LINK
      , nLinkType_(-1)
      #endif
      #ifdef SIMP_TETHER
      , hasTether_(-1)
      #endif
      , communicatorPtr_(&communicator)
   {
      setClassName("Simulation");
      Util::initStatic();
      Atom::initStatic();
      Analyzer::initStatic();

      if (!MPI::Is_initialized()) {
         UTIL_THROW("MPI not initialized on entry");
      }
      commitMpiTypes();

      // Set directory Id in FileMaster to MPI processor rank.
      fileMaster_.setDirectoryId(communicatorPtr_->Get_rank());

      // Set log file for processor n to a new file named "n/log"
      // Relies on initialization of FileMaster outputPrefix to "" (empty).
      fileMaster_.openOutputFile("log", logFile_);
      Log::setFile(logFile_);

      speciesManagerPtr_ = new SpeciesManager;
   }
   #endif

   /*
   * Constructor.
   */
   Simulation::Simulation()
    : iStep_(0), 
      nSystem_(1),
      speciesManagerPtr_(0),
      analyzerManagerPtr_(0),
      moleculeCapacity_(0),
      nAtomType_(-1),
      atomCapacity_(0)
      #ifndef SIMP_NOPAIR
      , maskedPairPolicy_(MaskBonded)
      #endif
      #ifdef SIMP_BOND
      , nBondType_(-1)
      , bondCapacity_(0)
      #endif
      #ifdef SIMP_ANGLE
      , nAngleType_(-1)
      , angleCapacity_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , nDihedralType_(-1)
      , dihedralCapacity_(0)
      #endif
      #ifdef SIMP_COULOMB
      , hasCoulomb_(-1)
      #endif
      #ifdef SIMP_EXTERNAL
      , hasExternal_(-1)
      #endif
      #ifdef MCMD_LINK
      , nLinkType_(-1)
      #endif
      #ifdef SIMP_TETHER
      , hasTether_(-1)
      #endif
      #ifdef UTIL_MPI
      , communicatorPtr_(0)
      #endif
   {
      setClassName("Simulation");
      Util::initStatic();
      Atom::initStatic();
      Analyzer::initStatic();
      speciesManagerPtr_ = new SpeciesManager;
   }

   /*
   * Destructor.
   */
   Simulation::~Simulation()
   {
      if (speciesManagerPtr_) {
         delete speciesManagerPtr_;
      }

      Atom::deallocate();

      #ifdef UTIL_MPI
      if (logFile_.is_open()) logFile_.close();
      #endif
   }

   #ifdef UTIL_MPI
   /*
   * Set an MPI job to read a single parameter file, from std::cin.
   */
   void Simulation::setIoCommunicator(MPI::Intracomm& communicator)
   {
      if (communicatorPtr_ == 0) {
         UTIL_THROW("No communicator was passed to constructor");
      } else 
      if (communicatorPtr_ != &communicator) {
         UTIL_THROW("ParamCcommunicator must be the one passed to constructor");
      }
      fileMaster_.setCommonControl();
      ParamComponent::setIoCommunicator(communicator);
   }

   /*
   * Set an MPI job to read a single parameter file, from std::cin.
   */
   void Simulation::setIoCommunicator()
   {  Simulation::setIoCommunicator(communicator()); }
   #endif

   /*
   * Read parameter file.
   */
   void Simulation::readParameters(std::istream& in)
   {
      // Preconditions
      assert(speciesManagerPtr_);

      readParamComposite(in, fileMaster_);

      read<int>(in, "nAtomType", nAtomType_);
      if (nAtomType_ <= 0) {
         UTIL_THROW("nAtomType must be > 0");
      }
      #ifdef SIMP_BOND
      nBondType_ = 0; // Default value
      readOptional<int>(in, "nBondType", nBondType_);
      if (nBondType_ < 0) {
         UTIL_THROW("nBondType must be >= 0");
      }
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_ = 0; // Default value
      readOptional<int>(in, "nAngleType", nAngleType_); 
      if (nAngleType_ < 0) {
         UTIL_THROW("nAngleType must be >= 0");
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_ = 0; // Default value 
      readOptional<int>(in, "nDihedralType", nDihedralType_); 
      if (nDihedralType_ < 0) {
         UTIL_THROW("nDihedralType must be >= 0");
      }
      #endif
      #ifdef SIMP_COULOMB
      hasCoulomb_ = 0;
      readOptional<int>(in, "hasCoulomb", hasCoulomb_); 
      if ((hasCoulomb_ != 0) && (hasCoulomb_ != 1)) {
         UTIL_THROW("hasCoulomb must be 0 or 1");
      }
      #endif
      #ifdef SIMP_EXTERNAL
      hasExternal_ = 0; // Default value 
      readOptional<int>(in, "hasExternal", hasExternal_); 
      if (hasExternal_ < 0) {
         UTIL_THROW("hasExternal must be >= 0");
      }
      #endif
      #ifdef MCMD_LINK
      nLinkType_ = 0; // Default value
      readOptional<int>(in, "nLinkType", nLinkType_); 
      if (nLinkType_ < 0) {
         UTIL_THROW("nLinkType must be >= 0");
      }
      #endif
      #ifdef SIMP_TETHER
      read<int>(in, "hasTether", hasTether_);
      if (hasTether_ < 0) {
         UTIL_THROW("hasTether must be >= 0");
      }
      #endif

      // Allocate and initialize array of AtomType objects
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
         #ifdef SIMP_COULOMB
         atomTypes_[i].setHasCharge(hasCoulomb_);
         #endif
      }

      // Read AtomType data from file
      readDArray<AtomType>(in, "atomTypes", atomTypes_, nAtomType_);

      #ifndef SIMP_NOPAIR
      read<MaskPolicy>(in, "maskedPairPolicy", maskedPairPolicy_);
      #endif

      readParamComposite(in, *speciesManagerPtr_);
      for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
         species(iSpecies).setId(iSpecies);
      }

      readParamComposite(in, random_);

      // Allocate and initialize all private arrays.
      initialize();

   }

   /*
   * Load internal state from an archive.
   */
   void Simulation::loadParameters(Serializable::IArchive &ar)
   {
      loadParamComposite(ar, fileMaster_);

      loadParameter<int>(ar, "nAtomType", nAtomType_);
      #ifdef SIMP_BOND
      loadParameter<int>(ar, "nBondType", nBondType_);
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_ = 0;
      loadParameter<int>(ar, "nAngleType", nAngleType_, false);
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_ = 0;
      loadParameter<int>(ar, "nDihedralType", nDihedralType_, false);
      #endif
      #ifdef SIMP_COULOMB
      hasCoulomb_ = 0;
      loadParameter<int>(ar, "hasCoulomb", hasCoulomb_, false);
      #endif
      #ifdef SIMP_EXTERNAL
      hasExternal_ = 0;
      loadParameter<int>(ar, "hasExternal", hasExternal_, false);
      #endif
      #ifdef MCMD_LINK
      nLinkType_ = 0;
      loadParameter<int>(ar, "nLinkType", nLinkType_, false);
      #endif
      #ifdef SIMP_TETHER
      hasTether_ = false;
      loadParameter<int>(ar, "hasTether", hasTether_, false);
      #endif

      // Allocate and load an array of AtomType objects
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      loadDArray<AtomType>(ar, "atomTypes", atomTypes_, nAtomType_);

      #ifndef SIMP_NOPAIR
      loadParameter<MaskPolicy>(ar, "maskedPairPolicy", maskedPairPolicy_);
      #endif
      loadParamComposite(ar, *speciesManagerPtr_);
      for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
         species(iSpecies).setId(iSpecies);
      }
      loadParamComposite(ar, random_);

      // Allocate and initialize all private arrays.
      initialize();
   }

   /*
   * Load internal state from an archive.
   */
   void Simulation::save(Serializable::OArchive &ar)
   {
      fileMaster_.save(ar);
      ar << nAtomType_;
      #ifdef SIMP_BOND
      ar << nBondType_;
      #endif
      #ifdef SIMP_ANGLE
      Parameter::saveOptional(ar, nAngleType_, (bool)nAngleType_);
      #endif
      #ifdef SIMP_DIHEDRAL
      Parameter::saveOptional(ar, nDihedralType_, (bool)nDihedralType_);
      #endif
      #ifdef SIMP_COULOMB
      Parameter::saveOptional(ar, hasCoulomb_, hasCoulomb_);
      #endif
      #ifdef SIMP_EXTERNAL
      Parameter::saveOptional(ar, hasExternal_, hasExternal_);
      #endif
      #ifdef MCMD_LINK
      Parameter::saveOptional(ar, nLinkType_, (bool)nLinkType_);
      #endif
      #ifdef SIMP_TETHER
      ar << hasTether_;
      Parameter::saveOptional(ar, hasTether_, hasTether_);
      #endif

      ar << atomTypes_;
      #ifndef SIMP_NOPAIR
      ar & maskedPairPolicy_;
      #endif
      (*speciesManagerPtr_).save(ar);
      random_.save(ar);
   }

   /*
   * Allocate and initialize all private data (private function).
   *
   * Allocates global arrays (molecules_, atoms_, bonds_, angles_) and the
   * arrays first<class>Ids_ of integers to species blocks. Initializes:
   *
   *   - Capacity values and first<class>Ptr_ addresses.
   *   - Integer ids for Species and Molecule objects.
   *   - Pointers between Species, Molecule, and Atom objects
   *   - Atom typeIds and all Bond and Angle objects.
   */
   void Simulation::initialize()
   {
      //Preconditions
      assert(nSpecies() > 0);
      if (nSpecies() <= 0) {
         UTIL_THROW("Error: nSpecies() <= 0 in Simulation::initialize()");
      }
      #ifdef SIMP_BOND
      if (nBondType_ < 0) {
         UTIL_THROW("Error: nBondType < 0 in Simulation::initialize()");
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_ < 0) {
         UTIL_THROW("Error: nAngleType < 0 in Simulation::initialize()");
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_ < 0) {
         UTIL_THROW("Error: nDihedralType < 0 in Simulation::initialize()");
      }
      #endif
      #ifdef MCMD_LINK
      if (nLinkType_ < 0) {
         UTIL_THROW("Error: nLinkType_ < 0 in Simulation::initialize()");
      }
      #endif

      Species *speciesPtr;
      int nAtom, iSpecies;
      int capacity;
      #ifdef SIMP_BOND
      int nBond;
      #endif
      #ifdef SIMP_ANGLE
      int nAngle;
      #endif
      #ifdef SIMP_DIHEDRAL
      int nDihedral;
      #endif

      // Allocate arrays dimensioned by nSpecies().
      reservoirs_.allocate(nSpecies());
      firstMoleculeIds_.allocate(nSpecies());
      firstAtomIds_.allocate(nSpecies());
      #ifdef SIMP_BOND
      if (nBondType_ > 0) {
         firstBondIds_.allocate(nSpecies());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_ > 0) {
         firstAngleIds_.allocate(nSpecies());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_ > 0) {
         firstDihedralIds_.allocate(nSpecies());
      }
      #endif

      // Count Molecules, Atoms and Groups.
      moleculeCapacity_ = 0;
      atomCapacity_ = 0;
      #ifdef SIMP_BOND
      bondCapacity_ = 0;
      #endif
      #ifdef SIMP_ANGLE
      angleCapacity_ = 0;
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralCapacity_ = 0;
      #endif
      for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
         speciesPtr = &species(iSpecies);

         // Check species id
         if (speciesPtr->id() != iSpecies) {
            UTIL_THROW("Inconsistent species ids");
         }
         //speciesPtr->setId(iSpecies);

         // Allocate reservoir for this species
         reservoirs_[iSpecies].allocate(speciesPtr->capacity());

         // Set indexes of first objects of the blocks for this species
         firstMoleculeIds_[iSpecies] = moleculeCapacity_;
         firstAtomIds_[iSpecies] = atomCapacity_;
         #ifdef SIMP_BOND
         if (nBondType_ > 0) {
            firstBondIds_[iSpecies] = bondCapacity_;
         }
         #endif
         #ifdef SIMP_ANGLE
         if (nAngleType_ > 0) {
            firstAngleIds_[iSpecies] = angleCapacity_;
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (nDihedralType_ > 0) {
            firstDihedralIds_[iSpecies] = dihedralCapacity_;
         }
         #endif

         // Increment total capacity values
         capacity = speciesPtr->capacity();
         nAtom = speciesPtr->nAtom();
         moleculeCapacity_ += capacity;
         atomCapacity_ += capacity*nAtom;
         #ifdef SIMP_BOND
         if (nBondType_ > 0) {
            nBond = speciesPtr->nBond();
            bondCapacity_ += capacity*nBond;
         }
         #endif
         #ifdef SIMP_ANGLE
         if (nAngleType_ > 0) {
            nAngle = speciesPtr->nAngle();
            angleCapacity_ += capacity*nAngle;
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (nDihedralType_ > 0) {
            nDihedral = speciesPtr->nDihedral();
            dihedralCapacity_ += capacity*nDihedral;
         }
         #endif
      }

      // Allocate global array of atoms (static member of Atom class).
      Atom::allocate(atomCapacity_, atoms_);

      // Allocate other global arrays (members of Simulation).
      molecules_.allocate(moleculeCapacity_);

      // Initialize all Atoms and Molecule objects.
      for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
         speciesPtr = &species(iSpecies);
         initializeSpecies(iSpecies);
      }

      #ifdef SIMP_BOND
      // Initialize bonds.
      if (nBondType_ > 0) {
         if (bondCapacity_ > 0) {
            bonds_.allocate(bondCapacity_);
         } else {
            bonds_.allocate(1);
         }
         for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
            initializeSpeciesBonds(iSpecies);
         }
      }
      #endif

      #ifdef SIMP_ANGLE
      // Initialize angles.
      if (nAngleType_ > 0) {
         if (angleCapacity_ > 0) {
            angles_.allocate(angleCapacity_);
         } else {
            angles_.allocate(1);
         }
         for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
            initializeSpeciesAngles(iSpecies);
         }
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Initialize dihedrals.
      if (nDihedralType_ > 0) {
         if (dihedralCapacity_ > 0) {
            dihedrals_.allocate(dihedralCapacity_);
         } else {
            dihedrals_.allocate(1);
         }
         for (iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {
            initializeSpeciesDihedrals(iSpecies);
         }
      }
      #endif

   }

   /*
   * Initialize all Molecule and Atom objects for one Species (private).
   *
   * This function creates associations between Species, Molecule, and
   * Atom objects for all molecules of one species, and sets atom typeIds.
   *
   * For each molecule, it sets the id, species pointer, nAtom, and the 
   * firstAtom pointer. The molecule id is only unique within each species.
   *
   * For each atom, it sets the molecule pointer and an integer typeId.
   *
   * This function also pushes all molecules of the species onto the
   * reservoir, pushing them in order of decreasing molecule id.
   */
   void Simulation::initializeSpecies(int iSpecies)
   {

      Species* speciesPtr;
      Molecule* moleculePtr;
      Atom* atomPtr;
      int iMol, iAtom;
      int capacity, nAtom;

      speciesPtr = &species(iSpecies);
      capacity   = speciesPtr->capacity();
      nAtom      = speciesPtr->nAtom();

      // Initialize pointers before loop
      moleculePtr = &molecules_[firstMoleculeIds_[iSpecies]];
      atomPtr     = &atoms_[firstAtomIds_[iSpecies]];

      // Loop over all molecules in Species
      for (iMol = 0; iMol < capacity; ++iMol) {

         // Initialize a Molecule
         moleculePtr->setId(iMol);
         moleculePtr->setSpecies(*speciesPtr);
         moleculePtr->setNAtom(nAtom);
         moleculePtr->setFirstAtom(*atomPtr);

         // Loop over atoms in a molecule, set molecule and atom TypeId
         for (iAtom = 0; iAtom < nAtom; ++iAtom) {
            atomPtr->setMolecule(*moleculePtr);
            atomPtr->setTypeId(speciesPtr->atomTypeId(iAtom));
            ++atomPtr;
         }

         ++moleculePtr;
      }

      // Push all molecules of this species onto the reservoir stack
      // Push on in reverse order, so that they pop off in sequence
      moleculePtr = &molecules_[firstMoleculeIds_[iSpecies] + capacity - 1];
      for (iMol = 0; iMol < capacity; ++iMol) {
         reservoirs_[iSpecies].push(*moleculePtr);
         --moleculePtr;
      }

   }

   #ifdef SIMP_BOND
   /*
   * Initialize all Bond objects for Molecules of one Species. (private)
   *
   * This functions assigns pointers to Atoms and bond types ids within a
   * contiguous block of Bond objects, and sets a pointer in each Molecule
   * to the first Bond in the associated block.
   */
   void Simulation::initializeSpeciesBonds(int iSpecies)
   {
      if (nBondType_ <= 0) {
         UTIL_THROW("nBondType_ must be positive");
      }

      Species* speciesPtr = 0;
      Molecule* moleculePtr = 0;
      Bond* bondPtr = 0;
      Atom* firstAtomPtr;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int iMol, iBond, atom0Id, atom1Id, type;
      int capacity, nBond;

      speciesPtr = &species(iSpecies);
      nBond = speciesPtr->nBond();
      capacity = speciesPtr->capacity();

      // Initialize pointers before loop
      moleculePtr = &molecules_[firstMoleculeIds_[iSpecies]];
      // bondPtr = &bonds_[firstBondIds_[iSpecies]];
      bondPtr = &bonds_[0] + firstBondIds_[iSpecies];

      // Loop over molecules in Species
      for (iMol = 0; iMol < capacity; ++iMol) {

         firstAtomPtr = &(moleculePtr->atom(0));
         moleculePtr->setFirstBond(*bondPtr);
         moleculePtr->setNBond(nBond);

         if (nBond > 0) {

            // Create bonds for a molecule
            for (iBond = 0; iBond < nBond; ++iBond) {

               // Get pointers to bonded atoms and bond type
               atom0Id  = speciesPtr->speciesBond(iBond).atomId(0);
               atom1Id  = speciesPtr->speciesBond(iBond).atomId(1);
               type     = speciesPtr->speciesBond(iBond).typeId();
               atom0Ptr = firstAtomPtr + atom0Id;
               atom1Ptr = firstAtomPtr + atom1Id;

               // Set fields of the Bond object
               bondPtr->setAtom(0, *atom0Ptr);
               bondPtr->setAtom(1, *atom1Ptr);
               bondPtr->setTypeId(type);

               #ifndef SIMP_NOPAIR
               // If MaskBonded, add each bonded atom to its partners Mask
               if (maskedPairPolicy_ == MaskBonded) {
                  atom0Ptr->mask().append(*atom1Ptr);
                  atom1Ptr->mask().append(*atom0Ptr);
               }
               #endif

               ++bondPtr;
            }

         }
         ++moleculePtr;
      }

   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Initialize all Angle objects for Molecules of one Species.
   *
   * This functions assigns pointers to Atoms and angle types ids within a
   * contiguous block of Angle objects, and sets a pointer in each Molecule
   * to the first Angle in the associated block.
   */
   void Simulation::initializeSpeciesAngles(int iSpecies)
   {

      if (nAngleType_ <= 0) {
         UTIL_THROW("nAngleType must be positive");
      }

      Species* speciesPtr = 0;
      Molecule* moleculePtr = 0;
      Angle* anglePtr = 0;
      Atom* firstAtomPtr, *atom0Ptr, *atom1Ptr, *atom2Ptr;
      int iMol, iAngle, atom0Id, atom1Id, atom2Id, type;
      int capacity, nAngle;

      speciesPtr = &species(iSpecies);
      capacity = speciesPtr->capacity();
      nAngle = speciesPtr->nAngle();

      // Initialize pointers before loop
      moleculePtr = &molecules_[firstMoleculeIds_[iSpecies]];
      // anglePtr = &angles_[firstAngleIds_[iSpecies]];
      anglePtr = &angles_[0] + firstAngleIds_[iSpecies];

      // Loop over molecules in Species
      for (iMol = 0; iMol < capacity; ++iMol) {

         firstAtomPtr = &(moleculePtr->atom(0));
         moleculePtr->setFirstAngle(*anglePtr);
         moleculePtr->setNAngle(nAngle);

         if (nAngle > 0) {

            // Create angles for a molecule
            for (iAngle = 0; iAngle < nAngle; ++iAngle) {

               // Get pointers to atoms spanning the angle and angle type
               atom0Id  = speciesPtr->speciesAngle(iAngle).atomId(0);
               atom1Id  = speciesPtr->speciesAngle(iAngle).atomId(1);
               atom2Id  = speciesPtr->speciesAngle(iAngle).atomId(2);
               type     = speciesPtr->speciesAngle(iAngle).typeId();
               atom0Ptr = firstAtomPtr + atom0Id;
               atom1Ptr = firstAtomPtr + atom1Id;
               atom2Ptr = firstAtomPtr + atom2Id;

               // Set fields of the Angle object
               anglePtr->setAtom(0, *atom0Ptr);
               anglePtr->setAtom(1, *atom1Ptr);
               anglePtr->setAtom(2, *atom2Ptr);
               anglePtr->setTypeId(type);

               ++anglePtr;

            }

         }

         ++moleculePtr;
      }

   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Initialize all Dihedral objects for Molecules of one Species.
   *
   * This functions assigns pointers to Atoms and Dihedral types ids within
   * a contiguous block of Dihedral objects, and sets a pointer in each 
   * Molecule to the first Dihedral in the associated block.
   */
   void Simulation::initializeSpeciesDihedrals(int iSpecies)
   {

      Species* speciesPtr = 0;
      Molecule* moleculePtr = 0;
      Dihedral* dihedralPtr = 0;
      Atom *firstAtomPtr, *atom0Ptr, *atom1Ptr, *atom2Ptr, *atom3Ptr;
      int iMol, iDihedral, atom0Id, atom1Id, atom2Id, atom3Id, type;
      int capacity, nDihedral;

      speciesPtr = &species(iSpecies);
      capacity   = speciesPtr->capacity();
      nDihedral  = speciesPtr->nDihedral();

      // Initialize pointers before loop
      moleculePtr = &molecules_[firstMoleculeIds_[iSpecies]];
      //dihedralPtr = &dihedrals_[firstDihedralIds_[iSpecies]];
      dihedralPtr = &dihedrals_[0] + firstDihedralIds_[iSpecies];

      // Loop over molecules in Species
      for (iMol = 0; iMol < capacity; ++iMol) {

         firstAtomPtr = &(moleculePtr->atom(0));
         moleculePtr->setFirstDihedral(*dihedralPtr);
         moleculePtr->setNDihedral(nDihedral);

         if (nDihedral > 0) {

            // Create dihedrals for a molecule
            for (iDihedral = 0; iDihedral < nDihedral; ++iDihedral) {

               // Get local indices for atoms and dihedral type
               atom0Id  = speciesPtr->speciesDihedral(iDihedral).atomId(0);
               atom1Id  = speciesPtr->speciesDihedral(iDihedral).atomId(1);
               atom2Id  = speciesPtr->speciesDihedral(iDihedral).atomId(2);
               atom3Id  = speciesPtr->speciesDihedral(iDihedral).atomId(3);
               type     = speciesPtr->speciesDihedral(iDihedral).typeId();

               // Calculate atom pointers
               atom0Ptr = firstAtomPtr + atom0Id;
               atom1Ptr = firstAtomPtr + atom1Id;
               atom2Ptr = firstAtomPtr + atom2Id;
               atom3Ptr = firstAtomPtr + atom3Id;

               // Set fields of the Dihedral object
               dihedralPtr->setAtom(0, *atom0Ptr);
               dihedralPtr->setAtom(1, *atom1Ptr);
               dihedralPtr->setAtom(2, *atom2Ptr);
               dihedralPtr->setAtom(3, *atom3Ptr);
               dihedralPtr->setTypeId(type);

               ++dihedralPtr;
            }
         }
         ++moleculePtr;
      }

   }
   #endif

   /*
   * Associate and allocate a System MoleculeSet for a specified Species.
   */
   void Simulation::allocateMoleculeSet(Util::ArraySet<Molecule>& set, 
                                        int speciesId) const
   {
      const Molecule* molecules = &molecules_[firstMoleculeIds_[speciesId]];
      int capacity = species(speciesId).capacity();
      set.allocate(molecules, capacity);
   }

   /*
   * Get a new molecule from a reservoir of unused Molecule objects.
   */ 
   Molecule& Simulation::getMolecule(int speciesId)
   {
      Molecule* ptr = &reservoirs_[speciesId].pop();  
      Activate::activate(*ptr); // activate all atoms in molecule
      return *ptr;
   }

   /*
   * Return a molecule to a reservoir of unused molecules.
   */ 
   void Simulation::returnMolecule(Molecule& molecule)
   {
      int speciesId = molecule.species().id();
      reservoirs_[speciesId].push(molecule);  
   }

   // Accessors

   /*
   * Get number of Species in this Simulation.
   */
   int Simulation::nSpecies() const
   {  return speciesManagerPtr_->size(); }

   /*
   * Get a reference to a Species object, identified by index.
   */
   Species& Simulation::species(int i)
   {  return (*speciesManagerPtr_)[i]; }

   /*
   * Get a const reference to a Species object, identified by index.
   */
   const Species& Simulation::species(int i) const
   {  return (*speciesManagerPtr_)[i]; }

   /*
   * Return true if valid, or throw an Exception.
   */
   bool Simulation::isValid() const
   {
      // Check indices of all Atom objects
      for (int iAtom = 0; iAtom < atomCapacity_; ++iAtom) {
         if (atoms_[iAtom].id() != iAtom) {
            UTIL_THROW("Error in atom Ids");
         }
      }

      // Check validity of all Species objects.
      for (int iSpecies=0; iSpecies < nSpecies(); ++iSpecies) {
         species(iSpecies).isValid();
      }

      // Declare and initialize pointers to atoms, bonds, etc.
      const Atom* atomPtr = &(atoms_[0]);
      #ifdef SIMP_BOND
      const Bond* bondPtr = 0;
      if (nBondType_ > 0) {
         bondPtr = &(bonds_[0]);
      }
      #endif
      #ifdef SIMP_ANGLE
      const Angle* anglePtr = 0;
      if (nAngleType_ > 0) {
         anglePtr = &(angles_[0]);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      const Dihedral* dihedralPtr = 0;
      if (nDihedralType_ > 0) {
          dihedralPtr  = &(dihedrals_[0]);
      }
      #endif
      const Molecule* moleculePtr = &(atomPtr->molecule());

      // Loop over molecular species
      for (int iSpecies = 0; iSpecies < nSpecies(); ++iSpecies) {

         // Declare / initialize variables for this Species
         const Species* speciesPtr = &(species(iSpecies));
         const int capacity = speciesPtr->capacity();
         const int nAtom = speciesPtr->nAtom();
         #ifdef SIMP_BOND
         int nBond = 0;
         if (nBondType_ > 0) {
            nBond = speciesPtr->nBond();
         }
         #endif
         #ifdef SIMP_ANGLE
         int nAngle = 0;
         if (nAngleType_ > 0) {
            nAngle = speciesPtr->nAngle();
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         int nDihedral = 0;
         if (nDihedralType_ > 0) {
            nDihedral = speciesPtr->nDihedral();
         }
         #endif

         // Loop over molecules allocated for a species
         for (int iMol = 0; iMol < capacity; ++iMol) {

            // Validate molecule
            if (&moleculePtr->species() != speciesPtr) {
               UTIL_THROW("Inconsistent molecule.species()");
            }
            if (moleculePtr->id() != iMol) {
               UTIL_THROW("Inconsistent molecule.id()");
            }

            // Validate atoms within a molecule
            {
               if (&moleculePtr->atom(0) != atomPtr) {
                  UTIL_THROW("Error in molecule::bond()");
               }

               int type;
               for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
   
                  // Validate pointers linking atom and molecule
                  if (&(atomPtr->molecule()) != moleculePtr) {
                     UTIL_THROW("Inconsistent atom.molecule()");
                  }
                  if (atomPtr != &(moleculePtr->atom(iAtom))) {
                     UTIL_THROW("Inconsistent molecule.atom()");
                  }
  
                  // Require that atom be active 
                  if (!atomPtr->isActive()) {
                     UTIL_THROW("Atom is inactive");
                  }

                  // Validate atom type, if species is not mutable
                  if (!speciesPtr->isMutable()) {
                     type = atomPtr->typeId();
                     if (type < 0 || type >= nAtomType_) {
                        UTIL_THROW("Invalid atom type id");
                     }
                     if (type != speciesPtr->atomTypeId(iAtom)) {
                        UTIL_THROW("Inconsistent atom type");
                     }
                  }


                  ++atomPtr;
               }
            }

            #ifdef SIMP_BOND
            if (nBondType_ > 0 && nBond > 0) {

               if (&moleculePtr->bond(0) != bondPtr) {
                  UTIL_THROW("Error in molecule::bond()");
               }
               if (moleculePtr->nBond() != nBond) {
                  UTIL_THROW("Inconsistent values of nBond");
               }

               const Atom* atom0Ptr;
               const Atom* atom1Ptr;
               int id0, id1, bondType;

               // Validate bonds within a molecule
               for (int iBond = 0; iBond < nBond; ++iBond) {

                  // Get data from associated SpeciesBond
                  id0 = speciesPtr->speciesBond(iBond).atomId(0);
                  id1 = speciesPtr->speciesBond(iBond).atomId(1);
                  bondType = speciesPtr->speciesBond(iBond).typeId();

                  // Validate SpeciesBond
                  if (id0 < 0 || id0 >= nAtom) {
                     UTIL_THROW("Invalid local atom id in a SpeciesBond");
                  }
                  if (id1 < 0 || id1 >= nAtom) {
                     UTIL_THROW("Invalid local atom id in a SpeciesBond");
                  }
                  if (bondType < 0 || bondType >= nBondType_) {
                     std::cout << "bondType   = " << bondType 
                              << std::endl;
                     std::cout << "nBondType_ = " << nBondType_ 
                              << std::endl;
                     UTIL_THROW("Invalid bondType in a SpeciesBond");
                  }

                  // Validate consistency of Bond and SpeciesBond
                  atom0Ptr = &(bondPtr->atom(0));
                  atom1Ptr = &(bondPtr->atom(1));
                  if (atom0Ptr != &(moleculePtr->atom(id0))) {
                     UTIL_THROW("Inconsistent atom handle for a Bond");
                  }
                  if (atom1Ptr != &(moleculePtr->atom(id1))) {
                     UTIL_THROW("Inconsistent atom handle for a Bond");
                  }
                  if (bondPtr->typeId() != bondType) {
                     UTIL_THROW("Inconsistent bond type index");
                  }

                  // Require that bond be active 
                  if (!bondPtr->isActive()) {
                     UTIL_THROW("Bond is inactive");
                  }

                  #ifndef SIMP_NOPAIR
                  // If MaskBonded, check that bonded atoms are masked
                  if (maskedPairPolicy_ == MaskBonded) {

                     if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {
                        UTIL_THROW("Missing masked partner");
                     }
                     if (!atom1Ptr->mask().isMasked(*atom0Ptr)) {
                        UTIL_THROW("Missing masked partner");
                     }

                  }
                  #endif

                  ++bondPtr;
               }
            }
            #endif

            #ifdef SIMP_ANGLE
            if (nAngleType_ > 0 && nAngle > 0) {

               if (&moleculePtr->angle(0) != anglePtr) {
                  UTIL_THROW("Error in molecule::angle()");
               }
               if (moleculePtr->nAngle() != nAngle) {
                  UTIL_THROW("Inconsistent values of nAngle");
               }

               const Atom* atom0Ptr = 0;
               const Atom* atom1Ptr = 0;
               const Atom* atom2Ptr = 0;
               int id0, id1, id2, angleType;

               // Validate angles within a molecule
               for (int iAngle = 0; iAngle < nAngle; ++iAngle) {

                  // Get data from associated SpeciesAngle
                  id0  = speciesPtr->speciesAngle(iAngle).atomId(0);
                  id1  = speciesPtr->speciesAngle(iAngle).atomId(1);
                  id2  = speciesPtr->speciesAngle(iAngle).atomId(2);
                  angleType = speciesPtr->speciesAngle(iAngle).typeId();

                  // Validate SpeciesAngle object
                  if (id0 < 0 || id0 >= nAtom) {
                     UTIL_THROW("Invalid local atom id in a SpeciesAngle");
                  }
                  if (id1 < 0 || id1 >= nAtom) {
                     UTIL_THROW("Invalid local atom id in a SpeciesAngle");
                  }
                  if (id2 < 0 || id2 >= nAtom) {
                     UTIL_THROW("Invalid local atom id in a SpeciesAngle");
                  }
                  if (angleType < 0 || angleType >= nAngleType_) {
                     UTIL_THROW("Invalid local atom id in a SpeciesAngle");
                  }

                  // Validate consistency of Angle and SpeciesAngle
                  atom0Ptr = &(anglePtr->atom(0));
                  atom1Ptr = &(anglePtr->atom(1));
                  atom2Ptr = &(anglePtr->atom(2));

                  if (atom0Ptr != &(moleculePtr->atom(id0))) {
                     UTIL_THROW("Inconsistent atom handle for an Angle");
                  }
                  if (atom1Ptr != &(moleculePtr->atom(id1))) {
                     UTIL_THROW("Inconsistent atom handle for an Angle");
                  }
                  if (atom2Ptr != &(moleculePtr->atom(id2))) {
                     UTIL_THROW("Inconsistent atom handle for an Angle");
                  }
                  if (anglePtr->typeId() != angleType) {
                     UTIL_THROW("Inconsistent angle type index");
                  }

                  // Require that angle be active 
                  if (!anglePtr->isActive()) {
                     UTIL_THROW("Angle is inactive");
                  }

                  ++anglePtr;
               }
            }
            #endif

            #ifdef SIMP_DIHEDRAL
            if (nDihedralType_ > 0 && nDihedral > 0) {

               const Atom* tAtomPtr;
               int tAtomId, dihedralType;

               if (&moleculePtr->dihedral(0) != dihedralPtr) {
                  UTIL_THROW("Error in molecule::dihedral()");
               }
               if (moleculePtr->nDihedral() != nDihedral) {
                  UTIL_THROW("Inconsistent values of nDihedral");
               }

               // Validate dihedrals within a molecule
               for (int iDihedral = 0; iDihedral < nDihedral; ++iDihedral) {

                  // Validate data from SpeciesDihedral
                  for (int tId = 0; tId < 4; ++tId) {
                     tAtomId =
                            speciesPtr->speciesDihedral(iDihedral).atomId(tId);
                     if (tId < 0 || tId >= nAtom)
                        UTIL_THROW("Invalid local atom id in SpeciesDihedral");

                     tAtomPtr = &(dihedralPtr->atom(tId));
                     if (tAtomPtr != &(moleculePtr->atom(tAtomId))) {
                        UTIL_THROW("Inconsistent atom handle for an Dihedral");
                     }
                  }
                  dihedralType = 
                               speciesPtr->speciesDihedral(iDihedral).typeId();
                  if (dihedralType < 0 || dihedralType >= nDihedralType_) {
                     UTIL_THROW("Invalid local atom id in a SpeciesDihedral");
                  }
                  if (dihedralType != dihedralPtr->typeId()) {
                     UTIL_THROW("Inconsistent dihedral type index");
                  }

                  // Require that dihedral be active 
                  if (!dihedralPtr->isActive()) {
                     UTIL_THROW("Dihedral is inactive");
                  }

                  ++dihedralPtr;
               }

            }
            #endif

            ++moleculePtr;
         }
      }

      return true;
   }

   // Mutators / Setters

   /*
   * Return the Species factory by reference.
   */
   Factory<Species>& Simulation::speciesFactory()
   {
      assert(speciesManagerPtr_);
      return speciesManagerPtr_->factory();
   }

   /*
   * Return the Analyzer factory by reference.
   */
   Factory<Analyzer>& Simulation::analyzerFactory()
   {
      assert(analyzerManagerPtr_);
      return analyzerManagerPtr_->factory();
   }

}
