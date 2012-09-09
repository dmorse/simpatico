#ifndef MCMD_SYSTEM_CPP
#define MCMD_SYSTEM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

// namespace McMd
#include "System.h"
#include "Simulation.h"
#include "serialize.h"

#include <mcMd/species/Species.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <mcMd/configIos/ConfigIo.h>
#include <mcMd/configIos/McConfigIo.h>
#include <mcMd/configIos/ConfigIoFactory.h>
#include <mcMd/trajectoryIos/TrajectoryIo.h>
#include <mcMd/trajectoryIos/TrajectoryIoFactory.h>

#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/PairFactory.h>
#endif
#include <mcMd/potentials/bond/BondFactory.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AngleFactory.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralFactory.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/potentials/link/LinkFactory.h>
#include <mcMd/links/LinkMaster.h>
#endif
#ifdef INTER_EXTERNAL
#include <mcMd/potentials/external/ExternalFactory.h>
#endif
#ifdef INTER_TETHER
#include <mcMd/potentials/tether/tetherFactory.h>
#include <mcMd/tethers/TetherMaster.h>
#endif

#ifdef MCMD_PERTURB
#include <mcMd/perturb/Perturbation.h>
#ifdef UTIL_MPI
#include <mcMd/perturb/ReplicaMove.h>
#endif
#endif

// namespace Util
#include <mcMd/util/FileMaster.h>
#include <util/param/Factory.h>
#include <util/archives/Serializable_includes.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   System::System()
    : moleculeSetsPtr_(0),
      boundaryPtr_(0),
      #ifdef MCMD_LINK
      linkMasterPtr_(0),
      #endif
      #ifdef INTER_TETHER
      tetherMasterPtr_(0),
      #endif
      simulationPtr_(0),
      energyEnsemblePtr_(0),
      boundaryEnsemblePtr_(0),
      #ifndef INTER_NOPAIR
      pairFactoryPtr_(0),
      #endif
      bondFactoryPtr_(0),
      #ifdef INTER_ANGLE
      angleFactoryPtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralFactoryPtr_(0),
      #endif
      #ifdef MCMD_LINK
      linkFactoryPtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalFactoryPtr_(0),
      #endif
      #ifdef INTER_TETHER
      tetherFactoryPtr_(0),
      #endif
      configIoPtr_(0),
      configIoFactoryPtr_(0),
      trajectoryIoFactoryPtr_(0),
      fileMasterPtr_(0),
      #ifdef MCMD_PERTURB
      perturbationPtr_(0),
      perturbationFactoryPtr_(0),
      #ifdef UTIL_MPI
      replicaMovePtr_(0),
      hasReplicaMove_(false),
      #endif
      #endif
      #ifndef INTER_NOPAIR
      pairStyle_(),
      #endif
      bondStyle_(),
      #ifdef INTER_ANGLE
      angleStyle_(),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStyle_(),
      #endif
      #ifdef MCMD_LINK
      linkStyle_(),
      #endif
      #ifdef INTER_EXTERNAL
      externalStyle_(),
      #endif
      #ifdef INTER_TETHER
      tetherStyle_(),
      #endif
      id_(0),
      isCopy_(false),
      createdFileMaster_(false)
      #ifdef MCMD_PERTURB
      , createdPerturbation_(false)
      , createdPerturbationFactory_(false)
      , createdReplicaMove_(false)
      #endif
   {
      moleculeSetsPtr_ = new DArray<MoleculeSet>;
      boundaryPtr_     = new Boundary;
      energyEnsemblePtr_   = new EnergyEnsemble;
      boundaryEnsemblePtr_ = new BoundaryEnsemble;
   }

   /*
   * Copy constructor.
   */
   System::System(const System& other)
    : ParamComposite(other),
      moleculeSetsPtr_(other.moleculeSetsPtr_),
      boundaryPtr_(other.boundaryPtr_),
      #ifdef MCMD_LINK
      linkMasterPtr_(other.linkMasterPtr_),
      #endif
      #ifdef INTER_TETHER
      tetherMasterPtr_(other.tetherMasterPtr_),
      #endif
      simulationPtr_(other.simulationPtr_),
      energyEnsemblePtr_(other.energyEnsemblePtr_),
      boundaryEnsemblePtr_(other.boundaryEnsemblePtr_),
      #ifndef INTER_NOPAIR
      pairFactoryPtr_(other.pairFactoryPtr_),
      #endif
      bondFactoryPtr_(other.bondFactoryPtr_),
      #ifdef INTER_ANGLE
      angleFactoryPtr_(other.angleFactoryPtr_),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralFactoryPtr_(other.dihedralFactoryPtr_),
      #endif
      #ifdef MCMD_LINK
      linkFactoryPtr_(other.linkFactoryPtr_),
      #endif
      #ifdef INTER_EXTERNAL
      externalFactoryPtr_(other.externalFactoryPtr_),
      #endif
      #ifdef INTER_TETHER
      tetherFactoryPtr_(other.tetherFactoryPtr_),
      #endif
      configIoPtr_(other.configIoPtr_),
      configIoFactoryPtr_(other.configIoFactoryPtr_),
      trajectoryIoFactoryPtr_(other.trajectoryIoFactoryPtr_),
      fileMasterPtr_(other.fileMasterPtr_),
      #ifdef MCMD_PERTURB
      perturbationPtr_(other.perturbationPtr_),
      perturbationFactoryPtr_(other.perturbationFactoryPtr_),
      #ifdef UTIL_MPI
      replicaMovePtr_(other.replicaMovePtr_),
      hasReplicaMove_(other.hasReplicaMove_),
      #endif
      #endif
      #ifndef INTER_NOPAIR
      pairStyle_(other.pairStyle_),
      #endif
      bondStyle_(other.bondStyle_),
      #ifdef INTER_ANGLE
      angleStyle_(other.angleStyle_),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStyle_(other.dihedralStyle_),
      #endif
      #ifdef MCMD_LINK
      linkStyle_(other.linkStyle_),
      #endif
      #ifdef INTER_EXTERNAL
      externalStyle_(other.externalStyle_),
      #endif
      #ifdef INTER_TETHER
      tetherStyle_(other.tetherStyle_),
      #endif
      id_(other.id_),
      isCopy_(true),
      createdFileMaster_(false)
      #ifdef MCMD_PERTURB
      , createdPerturbation_(false)
      , createdPerturbationFactory_(false)
      , createdReplicaMove_(false)
      #endif
   {}

   /*
   * Destructor.
   */
   System::~System()
   {
      if (!isCopy()) {
         if (moleculeSetsPtr_) {
            delete moleculeSetsPtr_;
         }
         if (boundaryPtr_) {
            delete boundaryPtr_;
         }
         #ifndef INTER_NOPAIR
         if (pairFactoryPtr_) {
            delete pairFactoryPtr_;
         }
         #endif
         if (bondFactoryPtr_) {
            delete bondFactoryPtr_;
         }
         #ifdef INTER_ANGLE
         if (angleFactoryPtr_) {
            delete angleFactoryPtr_;
         }
         #endif
         #ifdef INTER_DIHEDRAL
         if (dihedralFactoryPtr_) {
            delete dihedralFactoryPtr_;
         }
         #endif
         #ifdef MCMD_LINK
         if (linkFactoryPtr_) {
            delete linkFactoryPtr_;
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (externalFactoryPtr_) {
            delete externalFactoryPtr_;
         }
         #endif
         #ifdef MCMD_LINK
         if (linkMasterPtr_) {
            delete linkMasterPtr_;
         }
         #endif
         #ifdef INTER_TETHER
         if (tetherMasterPtr_) {
            delete tetherMasterPtr_;
         }
         if (tetherFactoryPtr_) {
            delete tetherFactoryPtr_;
         }
         #endif

         #ifdef MCMD_PERTURB
         if (perturbationPtr_ && createdPerturbationFactory_) {
            delete perturbationPtr_;
         }
         if (perturbationFactoryPtr_ && createdPerturbationFactory_) {
            delete perturbationFactoryPtr_;
         }
         #ifdef UTIL_MPI
         if (replicaMovePtr_ && createdReplicaMove_) {
            delete replicaMovePtr_;
         }
         #endif
         #endif
      
         if (energyEnsemblePtr_) {
            delete energyEnsemblePtr_;
         }
         if (boundaryEnsemblePtr_) {
            delete boundaryEnsemblePtr_;
         }
         if (configIoPtr_) {
            delete configIoPtr_;
         }
         if (configIoFactoryPtr_) {
            delete configIoFactoryPtr_;
         }
         if (trajectoryIoFactoryPtr_) {
            delete trajectoryIoFactoryPtr_;
         }
         if (fileMasterPtr_ && createdFileMaster_) {
            delete fileMasterPtr_;
         }
      }
   }

   /*
   * Read parameters from file.
   */
   void System::readParameters(std::istream &in)
   {

      // Only read parameters if this is not a copy.
      if (!isCopy()) {
         allocateMoleculeSets();
         readFileMaster(in);
         readPotentialStyles(in);
         #ifdef MCMD_LINK
         readLinkMaster(in);
         #endif
         #ifdef INTER_TETHER
         readTetherMaster(in);
         #endif
         readEnsembles(in);
      }

   }

   /**
   * If no FileMaster exists, create and initialize one. 
   */
   void System::readFileMaster(std::istream &in)
   {
      // Create FileMaster if necessary
      if (!fileMasterPtr_) {
         fileMasterPtr_ = new FileMaster();
         createdFileMaster_ = true;
         readParamComposite(in, *fileMasterPtr_);
      }
   }

   void System::readPotentialStyles(std::istream &in)
   {
      #ifndef INTER_NOPAIR
      read<std::string>(in, "pairStyle", pairStyle_);
      #endif

      if (simulation().nBondType() > 0) {
         read<std::string>(in, "bondStyle", bondStyle_);
      }

      #ifdef INTER_ANGLE
      if (simulation().nAngleType() > 0) {
         read<std::string>(in, "angleStyle", angleStyle_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (simulation().nDihedralType() > 0) {
         read<std::string>(in, "dihedralStyle", dihedralStyle_);
      }
      #endif

      #ifdef MCMD_LINK
      if (simulation().nLinkType() > 0) {
         read<std::string>(in, "linkStyle", linkStyle_);
      }
      #endif

      #ifdef INTER_EXTERNAL
      if (simulation().hasExternal()) {
         read<std::string>(in, "externalStyle", externalStyle_);
      }
      #endif

      #ifdef INTER_TETHER
      if (simulation().hasTether()) {
         read<std::string>(in, "tetherStyle", tetherStyle_);
      }
      #endif
   }

   /*
   * Create EnergyEnsemble and BoundaryEnsemble
   */
   void System::readEnsembles(std::istream &in)
   {
      readParamComposite(in, *energyEnsemblePtr_);
      readParamComposite(in, *boundaryEnsemblePtr_);
   }

   #ifdef MCMD_LINK
   void System::readLinkMaster(std::istream &in)
   {
      if (simulation().nLinkType() > 0) {
         linkMasterPtr_ = new LinkMaster();
         readParamComposite(in, *linkMasterPtr_);
      }
   }
   #endif 

   #ifdef INTER_TETHER
   void System::readTetherMaster(std::istream &in)
   {
      if (simulation().hasTether()) {
         tetherMasterPtr_ = new TetherMaster();
         readParamComposite(in, *tetherMasterPtr_);
      }
   }
   #endif 

   /*
   * Return a pointer to a new default ConfigIo.
   */
   ConfigIo* System::newDefaultConfigIo()
   {  return new McConfigIo(*this); }

   /*
   * Read configuration from a specific input stream.
   */
   void System::readConfig(std::istream &in)
   {
      if (!isEmpty()) removeAllMolecules();

      if (configIoPtr_ == 0) {
         configIoPtr_ = newDefaultConfigIo();
      }
      configIoPtr_->read(in);

      #ifdef UTIL_DEBUG
      isValid();
      #endif
   }

   /*
   * Write configuration to specified output stream.
   */
   void System::writeConfig(std::ostream &out)
   {
      if (configIoPtr_ == 0) {
         configIoPtr_ = newDefaultConfigIo();
      }
      configIoPtr_->write(out);
   }

   /* 
   * Save a System configuration to an archive.
   */
   void System::save(Serializable::OArchiveType& ar)
   {  ar & *this; }

   /* 
   * Load a System configuration from an archive.
   */
   void System::load(Serializable::IArchiveType& ar)
   {  ar & *this; }

   /*
   * Set System integer Id.
   *
   * Should be set immediately after the System is instantiated. If the System
   * is a member of the parent Simulation, setId() should be called in the
   * Simulation constructor.
   */
   void System::setId(int id)
   {  id_ = id; }

   /*
   * Set pointer to the parent Simulation.
   *
   * Should be set immediately after System instantiation. If this System
   * is a member of the parent Simulation, setSimulation(*this) should be
   * called in the Simulation constructor.
   */
   void System::setSimulation(Simulation &simulation)
   {
      assert(!simulationPtr_);
      simulationPtr_ = &simulation;
   }

   /*
   * Set pointer to a FileMaster.
   */
   void System::setFileMaster(FileMaster &fileMaster)
   {
      assert(!fileMasterPtr_);
      fileMasterPtr_ = &fileMaster;
   }

   // ConfigIoIo Management

   /*
   * Get the ConfigIo factory by reference.
   */
   Factory<ConfigIo>& System::configIoFactory()
   {
      if (!configIoFactoryPtr_) {
         configIoFactoryPtr_ = newDefaultConfigIoFactory();
      }
      return *configIoFactoryPtr_;
   }

   /*
   * Return a pointer to a new ConfigIoFactory.
   */
   Factory<ConfigIo>* System::newDefaultConfigIoFactory()
   {  return new ConfigIoFactory(*this); }

   /*
   * Set the ConfigIo, identified by subclass name.
   */
   void System::setConfigIo(std::string& classname)
   {
      if (!configIoFactoryPtr_) {
         configIoFactoryPtr_ = newDefaultConfigIoFactory();
      }
      ConfigIo* ptr = configIoFactoryPtr_->factory(classname);
      if (!ptr) {
         UTIL_THROW("Unrecognized ConfigIo subclass name");
      } 
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      configIoPtr_ = ptr;
   }

   // TrajectoryIo Management

   /*
   * Get the TrajectoryIo factory by reference.
   */
   Factory<TrajectoryIo>& System::trajectoryIoFactory()
   {
      if (!trajectoryIoFactoryPtr_) {
         trajectoryIoFactoryPtr_ = newDefaultTrajectoryIoFactory();
      }
      return *trajectoryIoFactoryPtr_;
   }

   /*
   * Return pointer to a new instance of default TrajectoryIo factory.
   */
   Factory<TrajectoryIo>* System::newDefaultTrajectoryIoFactory()
   {  return new TrajectoryIoFactory(*this); }

   #ifdef MCMD_PERTURB
   /*
   * Set this system to expect a Perturbation in the param file.
   * Also creates the perturbation factory.
   */
   void System::setExpectPerturbation()
   {
      // create perturbation factory
      perturbationFactoryPtr_ = newDefaultPerturbationFactory();
      createdPerturbationFactory_ = true;

      if (hasPerturbation()) {
         UTIL_THROW("A Perturbation is already set");
      }
      expectPerturbationParam_ = true;
   }

   void System::readPerturbation(std::istream& in) 
   {

      // Create Perturbation and read object, if required.
      if (!hasPerturbation() and expectPerturbationParam_) {
         std::string className;
         bool        isEnd;
         perturbationPtr_ = 
            perturbationFactoryPtr_->readObject(in, *this, className, isEnd);
         if (!perturbationPtr_) {
            std::string msg = "Unrecognized Perturbation subclass name ";
            msg += className;
            UTIL_THROW(msg.c_str());
         }
         createdPerturbation_ = true;
      }
   }
   #ifdef UTIL_MPI  
   void System::readReplicaMove(std::istream& in) 
   {
      if (hasPerturbation()) {
          read<bool>(in, "hasReplicaMove", hasReplicaMove_);
          if (hasReplicaMove_) {
             replicaMovePtr_ = new ReplicaMove(*this);
             readParamComposite(in, *replicaMovePtr_);
             // replicaMovePtr_->initialize();
          }
          createdReplicaMove_ = true;
      } else {
         hasReplicaMove_ = false;
      }
   }
   #endif // UTIL_MPI
   #endif // MCMD_PERTURB

   /*
   * Allocate and initialize a MoleculeSet for each Species.
   */
   void System::allocateMoleculeSets()
   {
      // Preconditions
      assert(simulationPtr_);

      // Allocate an array of nSpecies empty MoleculeSet objects.
      int nSpecies = simulationPtr_->nSpecies();
      moleculeSetsPtr_->allocate(nSpecies);

      // Allocate and initialize the MoleculeSet for each species.
      for (int i=0; i < nSpecies; ++i) {
         simulation().allocateMoleculeSet((*moleculeSetsPtr_)[i], i);
      }
   }

   /*
   * Add a molecule to this System.
   */
   void System::addMolecule(Molecule& molecule)
   {
      assert(moleculeSetsPtr_);

      int speciesId = molecule.species().id();
      if ( speciesId >= moleculeSetsPtr_->capacity()) {
         UTIL_THROW("Error: speciesId >= moleculeSetsPtr_->capacity()");
      }
      (*moleculeSetsPtr_)[speciesId].append(molecule);
      molecule.setSystem(*this);

      // notify observers
      notifyMoleculeSetObservers();
   }

   /*
   * Remove a molecule from this System.
   */
   void System::removeMolecule(Molecule& molecule)
   {
      assert(moleculeSetsPtr_);

      if ( &molecule.system() != this) {
         UTIL_THROW("Attempt to remove Molecule that is not in this System.");
      }
      int speciesId = molecule.species().id();
      (*moleculeSetsPtr_)[speciesId].remove(molecule);
      molecule.unsetSystem();

      // notify observers
      notifyMoleculeSetObservers();
   }

   /*
   * Is this System empty (i.e., devoid of Molecules) ?
   */
   bool System::isEmpty() const
   {
      for (int i = 0; i < simulation().nSpecies(); ++i) {
         if (nMolecule(i) != 0) return false;
      }
      return true;
   }

   /*
   * Return all molecules of all Species to their reservoirs.
   *
   * Ifdef MCMD_LINK, it also clears the LinkMaster.
   */
   void System::removeAllMolecules()
   {
      Species*  speciesPtr;
      Molecule* molPtr;
      int       iSpecies, nMol;
      for (iSpecies = 0; iSpecies < simulation().nSpecies(); ++iSpecies) {
         speciesPtr = &simulation().species(iSpecies);
         nMol  = nMolecule(iSpecies);
         while (nMol > 0) {
            molPtr = &molecule(iSpecies, nMol - 1);
            removeMolecule(*molPtr);
            speciesPtr->reservoir().push(*molPtr);
            nMol = nMolecule(iSpecies);
         }
      }

      #ifdef MCMD_LINK
      linkMaster().clear();
      #endif

      // notify observers
      notifyMoleculeSetObservers();
   }

   /*
   * Get the index of a molecule within its Species in this System.
   */
   int System::moleculeId(const Molecule& molecule) const
   {
      assert(moleculeSetsPtr_);

      int speciesId = molecule.species().id();
      return (*moleculeSetsPtr_)[speciesId].index(molecule);
   }

   /* 
   * Get a randomly chosen molecule of a specified Species.
   */
   Molecule& System::randomMolecule(int speciesId)
   {
      int nMol, moleculeId;
      nMol  = nMolecule(speciesId); 
      if (nMol <= 0) {
         UTIL_THROW("Number of molecules in species <= 0");
      }
      moleculeId = simulation().random().uniformInt(0, nMol);
      return molecule(speciesId, moleculeId); 
   }

   /*
   * Return the total number of atoms in this System.
   */
   int System::nAtom() const
   {
      int sum = 0;
      for (int i = 0; i < simulation().nSpecies(); ++i) {
         sum += nMolecule(i)*simulation().species(i).nAtom();
      }
      return sum;
   }

   #ifndef INTER_NOPAIR
   /*
   * Return the PairFactory by reference.
   */
   PairFactory& System::pairFactory()
   {
      if (!pairFactoryPtr_) {
         pairFactoryPtr_ = new PairFactory;
      }
      assert(pairFactoryPtr_);
      return *pairFactoryPtr_;
   }

   /*
   * Get the pair style string.
   */
   std::string System::pairStyle() const
   {  return pairStyle_;  }
   #endif

   /*
   * Return the BondFactory by reference.
   */
   Factory<BondPotential>& System::bondFactory()
   {
      if (!bondFactoryPtr_) {
         bondFactoryPtr_ = new BondFactory(*this);
      }
      assert(bondFactoryPtr_);
      return *bondFactoryPtr_;
   }

   /*
   * Get the bond style string.
   */
   std::string System::bondStyle() const
   {  return bondStyle_;  }

   #ifdef INTER_ANGLE
   /*
   * Return the AngleFactory by reference.
   */
   Factory<AnglePotential>& System::angleFactory()
   {
      if (angleFactoryPtr_ == 0) {
         angleFactoryPtr_ = new AngleFactory(*this);
      }
      assert(angleFactoryPtr_);
      return *angleFactoryPtr_;
   }

   /*
   * Get the angle style string.
   */
   std::string System::angleStyle() const
   {  return angleStyle_;  }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Return the DihedralFactory by reference.
   */
   Factory<DihedralPotential>& System::dihedralFactory()
   {
      if (dihedralFactoryPtr_ == 0) {
         dihedralFactoryPtr_ = new DihedralFactory(*this);
      }
      assert(dihedralFactoryPtr_);
      return *dihedralFactoryPtr_;
   }

   /*
   * Get the dihedral style string.
   */
   std::string System::dihedralStyle() const
   {  return dihedralStyle_;  }

   #endif

   #ifdef MCMD_LINK
   /*
   * Return the Link factory by reference.
   */
   Factory<BondPotential>& System::linkFactory()
   {
      if (linkFactoryPtr_ == 0) {
         linkFactoryPtr_ = new LinkFactory(*this);
      }
      assert(linkFactoryPtr_);
      return *linkFactoryPtr_;
   }

   /*
   * Get the link style string.
   */
   std::string System::linkStyle() const
   {  return linkStyle_;  }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Return the ExternalFactory by reference.
   */
   Factory<ExternalPotential>& System::externalFactory()
   {
      if (externalFactoryPtr_ == 0) {
         externalFactoryPtr_ = new ExternalFactory(*this);
      }
      assert(externalFactoryPtr_);
      return *externalFactoryPtr_;
   }

   /*
   * Get the external style string.
   */
   std::string System::externalStyle() const
   {  return externalStyle_;  }
   #endif

   #ifdef INTER_TETHER
   /*
   * Return the TetherFactory by reference.
   */
   TetherFactory& System::tetherFactory()
   {
      if (tetherFactoryPtr_ == 0) {
         tetherFactoryPtr_ = new TetherFactory;
      }
      assert(tetherFactoryPtr_);
      return *tetherFactoryPtr_;
   }

   /*
   * Get the tether style string.
   */
   std::string System::tetherStyle() const
   {  return tetherStyle_;  }
   #endif

   /*
   * Check validity of all data structures.
   */
   bool System::isValid() const
   {
      Molecule* molPtr;
      int       iSpecies, iMolecule;
      if (!simulationPtr_) {
         UTIL_THROW("Null simulationPtr_");
      }
      if (!moleculeSetsPtr_) {
         UTIL_THROW("Null moleculeSetsPtr_");
      }
      for (iSpecies =0; iSpecies < simulationPtr_->nSpecies(); ++iSpecies) {
         (*moleculeSetsPtr_)[iSpecies].isValid();
         for (iMolecule = 0; iMolecule < nMolecule(iSpecies); ++iMolecule) {
            molPtr = &(*moleculeSetsPtr_)[iSpecies][iMolecule]; 
            if (&molPtr->system() != this) {
               UTIL_THROW("Error in pointer to System in a Molecule");
            }
         }
      }

      #ifdef MCMD_LINK
      if (simulation().nLinkType() > 0) {
         if (!linkMasterPtr_) {
            UTIL_THROW("Null linkMasterPtr_");
         }
         linkMasterPtr_->isValid();
      }
      #endif

      #ifdef INTER_TETHER
      if (simulation().hasTether() > 0) {
         if (!tetherMasterPtr_) {
            UTIL_THROW("Null tetherMasterPtr_");
         }
         tetherMasterPtr_->isValid();
      }
      #endif

      return true;
   }

}
#endif
