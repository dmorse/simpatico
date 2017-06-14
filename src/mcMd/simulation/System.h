#ifndef MCMD_SYSTEM_H
#define MCMD_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <util/param/ParamComposite.h>        // base class

#include <util/boundary/Boundary.h>           // member (typedef)
#include <mcMd/chemistry/Molecule.h>          // member template parameter

#include <util/containers/DArray.h>           // member template
#include <util/containers/ArraySet.h>         // member template
#include <util/containers/PArrayIterator.h>   // inline function begin()

#include <iostream>
#include <string>
#include <set>

class SystemTest;

namespace Util 
{ 
   template <typename T> class Factory;
   class EnergyEnsemble;
   class BoundaryEnsemble;
   class FileMaster;
}

namespace McMd
{

   using namespace Util;

   /**
    * Observer interface. Classes that need to be notified
    * upon addition/removal of molecules can derive from this class
    * and register using System.addMoleculeSetObserver()
    */
   class MoleculeSetObserver 
   {
   public:

      virtual ~MoleculeSetObserver(){};
   
      virtual void notifyMoleculeSetChanged() = 0;
   };

   class Simulation;
   class ConfigIo;
   class TrajectoryReader;
   class PairFactory;
   #ifdef SIMP_BOND
   class BondPotential;
   #endif
   #ifdef SIMP_ANGLE
   class AnglePotential;
   #endif
   #ifdef SIMP_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef SIMP_COULOMB
   class CoulombFactory;
   #endif
   #ifdef SIMP_EXTERNAL
   class ExternalPotential;
   #endif
   #ifdef MCMD_LINK
   class LinkPotential;
   class LinkMaster;
   #endif 
   #ifdef SIMP_TETHER
   class TetherFactory;
   class TetherMaster;
   #endif 
   #ifdef MCMD_PERTURB
   class Perturbation;
   #ifdef UTIL_MPI
   class ReplicaMove;
   #endif
   #endif
   
   /**
   * A set of interacting Molecules enclosed by a Boundary.
   *
   * A System has:
   *
   *  - a Boundary, which defines the System dimensions.
   *  - an ArraySet of Molecule objects of each Species.
   *  - an EnergyEnsemble and a BoundaryEnsemble
   *  - a ConfigIo, which can read and write a configuration file
   *  - a FileMaster, to manage associated input and output files
   *
   * A System must be associated with a parent Simulation, which owns
   * many of the data structures used by a System. 
   *
   * The MdSystem and McSystem subclasses of System are designed for 
   * use in MD and MC simulations, respectively, and provide methods 
   * to evaluate energies and forces.
   *
   * \ingroup McMd_System_Module
   */
   class System : public ParamComposite
   {
   
   public:

      // Typedefs

      /// A set of molecules of one Species in a System.
      typedef ArraySet<Molecule>             MoleculeSet;

      /// Iterator for a MoleculeSet.
      typedef PArrayIterator<Molecule>       MoleculeIterator;

      /// Const Iterator for a MoleculeSet.
      typedef ConstPArrayIterator<Molecule>  ConstMoleculeIterator;

      // Methods

      /// Default constructor. 
      System();
 
      /// Copy constructor. 
      System(const System& other);
 
      /// Destructor.   
      virtual ~System();
 
      /// \name Initialization
      //@{
    
      /** 
      * Set the integer Id for this System.
      *
      * Set id to zero for Simulation objects with one System.
      *
      * \param Id integer System index.
      */
      void setId(int Id);
  
      /** 
      * Set the parent Simulation. 
      * 
      * \param simulation parent Simulation object.
      */
      void setSimulation(Simulation& simulation);

      /** 
      * Set the FileMaster. 
      * 
      * Is normally used to set the FileMaster to that of the parent 
      * Simulation.
      * 
      * \param filemaster FileMaster object.
      */
      void setFileMaster(FileMaster& filemaster);

      /**
      * Read parameter file.
      *
      * Only reads parameters if this System is not a copy (i.e.,
      * was not constructed with the copy constructor). If it is a 
      * copy, this function does nothing and returns normally.
      *
      * \param in pararameter file input stream
      */
      virtual void readParameters(std::istream& in);
 
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state from an archive.
      *
      * \param ar output/saving archive
      */
      void saveParameters(Serializable::OArchive &ar);

      //@}
      /// \name Config File IO
      //@{

      /**
      * Get the configuration file reader/writer factory by reference.
      */
      Factory<ConfigIo>& configIoFactory();

      /**
      * Create a new configuration file reader/writer.
      *
      * This function allows one to choose from among several subclasses
      * of ConfigIo, identified by subclass name. The implementation
      * uses a Factory<ConfigIo> object to instantiate a new object. 
      * If setConfigIoFactory() has not been called, an instance of
      * the default class ConfigIoFactory is created and used.
      *
      * \param classname name of desired ConfigIo subclass.
      */
      void setConfigIo(std::string& classname);

      /**
      * Read system configuration from file.
      *
      * The configuration file contains the dimensions of the
      * Boundary, the number of molecules of each Species, and 
      * the atomic positions for all atoms of all molecules in 
      * this System. 
      *
      * This function uses a ConfigIo object that is registered
      * with this System. If a ConfigIo object has not been
      * registered by calling setConfigIo(std::string&), this 
      * function creates and uses an instance of the ConfigIo 
      * base class.
      *
      * Precondition: The System and its Parent Simulation must
      * have been initialized by calling their readParam function. 
      *
      * \param in configuration file input stream
      */
      virtual void readConfig(std::istream& in);

      /**
      * Write system configuration to a specified ostream.
      *
      * Like readConfig(), this function will create and use 
      * a ConfigIo object if none has been registered 
      * previously. 
      *
      * \param out configuration file output stream
      */
      void writeConfig(std::ostream& out);

      /**
      * Load configuration.
      *
      * \param ar input/loading archive
      */
      virtual void loadConfig(Serializable::IArchive& ar);

      /**
      * Save configuration.
      *
      * \param ar output/save archive
      */
      void saveConfig(Serializable::OArchive& ar);

      //@}
      /// \name Trajectory File IO
      //@{

      /**
      * Get the trajectory reader/writer factory by reference.
      */
      Factory<TrajectoryReader>& trajectoryReaderFactory();

      //@}
      /// \name Molecule Set Mutators
      //@{
      
      /**
      * Add a Molecule to this System.  
      *
      * This function adds a Molecule to the set of Molecules of the same
      * Species in this System, and calls Molecule::setSystem(*this). 
      * 
      * The molecule to be added should first be popped off the Species 
      * reservoir or removed from another System.
      *
      * \param molecule Molecule to be added to this System.
      */
      void addMolecule(Molecule& molecule);

      /**
      * Remove a specific molecule from this System.  
      *
      * This function removes a Molecule from the set of molecules of the
      * same Species in this System, and calls Molecule::unsetSystem().
      *
      * The removed Molecule should be pushed onto the Species reservoir
      * or added to another System.
      *
      * \param molecule Molecule to be removed from this System.
      */
      void removeMolecule(Molecule& molecule); 

      /**
      * Remove all molecules from this System.
      *
      * Remove all molecules of every species from this System, and 
      * push each onto the reservoir for the appropriate Species.
      */
      void removeAllMolecules();

      /**
      * Subscribe to moleculeSet change signal
      *
      * \param observer the observer
      */
      void subscribeMoleculeSetChange(MoleculeSetObserver& observer);

      /**
      * Unsubscribe from moleculeSet change signal
      *
      * \param observer the observer
      */
      void unsubscribeMoleculeSetChange(MoleculeSetObserver& observer);

      //@}
      /// \name Molecule Set Accessors
      //@{
      
      /** 
      * Get the number of molecules of one Species in this System.
      *
      * \param speciesId integer Id for a Species.
      * \return number of molecules of specified Species in this System.
      */
      int nMolecule(int speciesId) const;

      /**
      * Return the total number of atoms in this System.
      */
      int nAtom() const;

      /**
      * Is this an empty System (i.e., one with no molecules) ?
      */
      bool isEmpty() const;

      /** 
      * Get the index of a Molecule within its Species in this System.
      *
      * This function returns the current (mutable) index of the molecule 
      * within the set of molecules of the same Species in this System, 
      * in the range 0 <= moleculeId < nMolecule(speciesId). This is 
      * the index required as the second argument of molecule(int, int).
      *
      * Note: The id returned by this function is not the same as the
      * id returned by Molecule::id(), which is a permanent identifier
      * for each molecule that is set immediately after allocation.
      *
      * \param molecule Molecule object of interest.
      * \return index for molecule within its Species and System.
      */
      int moleculeId(const Molecule& molecule) const;

      /** 
      * Get a specific Molecule in this System, by integer index.
      *
      * The moleculeId must be in range 0 <= moleculeId < nMolecule(speciesId).
      * The index associated with a molecule is mutable, and can change when
      * another molecule of the same Species is removed from this System. The
      * value of this index is returned by System::moleculeId().
      * 
      * \param speciesId  integer id of the desired Species
      * \param moleculeId integer id of molecule within Species and System
      * \return reference to the molecule
      */
      Molecule& molecule(int speciesId, int moleculeId);

      /** 
      * Get a random Molecule of a specified species in this System.
      *
      * \param speciesId  integer id of the desired Species.
      * \return const reference to the chosen molecule
      */
      Molecule& randomMolecule(int speciesId);

      /**
      * Initialize an iterator for molecules of one species in this System.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, MoleculeIterator& iterator);

      /**
      * Initialize a const iterator for molecules of one species in this System.
      *
      * \param speciesId integer Id for the desired Species (input)
      * \param iterator  molecule iterator (output)
      */
      void begin(int speciesId, ConstMoleculeIterator& iterator) const;

      //@}
      /// \name Potential Energy Factories and Styles
      //@{
    
      #ifndef SIMP_NOPAIR 
      /**
      * Get the PairFactory by reference.
      */
      PairFactory& pairFactory();

      /**
      * Return nonbonded pair style string.
      */
      std::string pairStyle() const;
      #endif

      #ifdef SIMP_BOND
      /**
      * Get the associated Factory<BondPotential> by reference.
      */
      Factory<BondPotential>& bondFactory();

      /**
      * Return covalent bond style string.
      */
      std::string bondStyle() const;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get the associated AngleFactory by reference.
      */
      Factory<AnglePotential>& angleFactory();

      /**
      * Return angle potential style string.
      */
      std::string angleStyle() const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get the associated Dihedral Factory by reference.
      */
      Factory<DihedralPotential>& dihedralFactory();

      /**
      * Return dihedral potential style string.
      */
      std::string dihedralStyle() const;
      #endif

      #ifdef SIMP_COULOMB
      /**
      * Get the associated Coulomb Factory by reference.
      */
      CoulombFactory& coulombFactory();

      /**
      * Return coulomb potential style string.
      */
      std::string coulombStyle() const;
      #endif

      #ifdef SIMP_EXTERNAL
      /**
      * Get the associated ExternalPotential factory by reference.
      */
      Factory<ExternalPotential>& externalFactory();

      /**
      * Return external potential style string.
      */
      std::string externalStyle() const;
      #endif

      #ifdef MCMD_LINK
      /**
      * Get the associated Link factory by reference.
      */
      Factory<BondPotential>& linkFactory();

      /**
      * Return link potential style string.
      */
      std::string linkStyle() const;

      /**
      * Get the LinkMaster by reference.
      */
      LinkMaster& linkMaster() const;
      #endif

      #ifdef SIMP_TETHER
      /**
      * Get the associated TetherFactory by reference.
      */
      TetherFactory& tetherFactory();

      /**
      * Return tether potential style string.
      */
      std::string tetherStyle() const;

      /**
      * Get the TetherMaster by reference.
      */
      TetherMaster& tetherMaster() const;
      #endif

      //@}
      #ifdef MCMD_PERTURB
      /// \name Free Energy Perturbation Theory
      //@{

      /**
      * Get the perturbation factory by reference.
      */
      Factory<Perturbation>& perturbationFactory();

      /**
      * Set to expect a Perturbation in the parameter file.
      */
      void setExpectPerturbation();

      /**
      * Return true if we expect a perturbation
      */
      bool expectPerturbation() const;

      /**
      * Does this system have an associated Perturbation?
      */
      bool hasPerturbation() const;

      /**
      * Get the associated Perturbation by reference.
      */
      Perturbation& perturbation() const;
     
      #ifdef UTIL_MPI 
      /**
      * Does this system have an associated ReplicaMove?
      */
      bool hasReplicaMove() const;
      
      /**
      * Get the associated ReplicaMove by reference.
      */
      ReplicaMove& replicaMove() const;
      #endif // UTIL_MPI

      //@}
      #endif // MCMD_PERTURB
      /// \name Accessors (Miscellaneous)
      //@{
    
      /**
      * Get integer index for this System.
      */ 
      int id() const;
  
      /**
      * Get the parent Simulation by reference.
      */ 
      Simulation& simulation() const;

      /**
      * Get the Boundary by reference.
      */ 
      Boundary& boundary() const;

      /**
      * Get the EnergyEnsemble by reference.
      */ 
      EnergyEnsemble& energyEnsemble() const;

      /**
      * Get the BoundaryEnsemble by reference.
      */ 
      BoundaryEnsemble& boundaryEnsemble() const;

      /**
      * Get the associated FileMaster by reference.
      */ 
      FileMaster& fileMaster() const;

      /**
      * Was this System instantiated with the copy constructor?
      */
      bool isCopy() const;

      /**
      * Return true if valid, or throw Exception. 
      */
      virtual bool isValid() const;

      //@}

   protected:

      /**
      * Get the maximum Boundary by reference.
      */
      Boundary& maxBoundary() const;

      /**
      * Return a pointer to a new default ConfigIo.
      */
      virtual ConfigIo* newDefaultConfigIo();

      /**
      * Return a pointer to a new default ConfigIoFactory.
      */
      virtual Factory<ConfigIo>* newDefaultConfigIoFactory();

      /**
      * Return a pointer to a new default TrajectoryReaderFactory.
      */
      virtual Factory<TrajectoryReader>* newDefaultTrajectoryReaderFactory();

      #ifdef MCMD_PERTURB
      /**
      * Return a pointer to the default perturbation Factory.
      */
      virtual Factory<Perturbation>* newDefaultPerturbationFactory()
      {  return 0; }

      /**
      * Read the perturbation parameter block (if any)
      *
      * \param in input parameter stream
      */
      void readPerturbation(std::istream& in);

      /**
      * Load the perturbation parameter block (if any)
      *
      * \param ar input/saving archive
      */
      void loadPerturbation(Serializable::IArchive& ar);

      /**
      * Save the perturbation parameter block (if any)
      *
      * \param ar output/saving archive
      */
      void savePerturbation(Serializable::OArchive& ar);

      #ifdef UTIL_MPI
      /**
      * Read the ReplicaMove parameter block (if any)
      *
      * \param in input parameter stream
      */
      void readReplicaMove(std::istream& in);

      /**
      * Read the ReplicaMove parameter block (if any)
      *
      * \param ar input/loading archive
      */
      void loadReplicaMove(Serializable::IArchive& ar);

      /**
      * Save the ReplicaMove parameter block (if any)
      *
      * \param ar output/saving archive
      */
      void saveReplicaMove(Serializable::OArchive& ar);
      #endif // UTIL_MPI
      #endif // MCMD_PERTURB

      /**
      * Allocate and initialize molecule sets for all species.
      *
      * This function is called within the Simulation::initialize() 
      * private member function to allocate and initialize an 
      * array of MoleculeSet objects for all Species for this 
      * System.
      *
      * Preconditions: This System must be associated with a parent
      * Simulation, and all Species objects must have been initialized
      * by calling SpeciesManager::readParameters().
      */
      void allocateMoleculeSets();

      /**
      * Read FileMaster parameters, if none yet exists.
      *
      * If no FileMaster exists, this function creates one and 
      * reads paramters to initialize one. If there is already
      * a FileMaster, it does nothing.
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readFileMaster(std::istream& in);

      /**
      * Load FileMaster data from archive, if necessary.
      *
      * \param ar input/loading archive
      */
      void loadFileMaster(Serializable::IArchive& ar);

      /**
      * If necessary, save FileMaster to archive.
      *
      * \param ar output/saving archive
      */
      void saveFileMaster(Serializable::OArchive& ar);

      /**
      * Read potential style parameter strings.
      *
      * \param in input parameter stream
      */
      void readPotentialStyles(std::istream& in);

      /**
      * Load potential style strings from an archive.
      *
      * \param ar input/loading archive
      */
      void loadPotentialStyles(Serializable::IArchive& ar);

      /**
      * Save potential style strings.
      *
      * \param ar output/saving archive
      */
      void savePotentialStyles(Serializable::OArchive& ar);

      /**
      * Read energy and boundary ensemble parameters.
      *
      * \param in input parameter stream
      */
      void readEnsembles(std::istream& in);

      /**
      * Load energy and boundary ensembles from archive.
      *
      * \param ar input/loading archive
      */
      void loadEnsembles(Serializable::IArchive& ar);

      /**
      * Save energy and boundary ensembles.
      *
      * \param ar output/saving archive
      */
      void saveEnsembles(Serializable::OArchive& ar);

      #ifdef MCMD_LINK
      /**
      * Read the LinkMaster parameters.
      *
      * \param in input parameter stream
      */
      void readLinkMaster(std::istream& in);

      /**
      * Load the LinkMaster.
      *
      * \param ar input archive.
      */
      void loadLinkMaster(Serializable::IArchive& ar);

      /**
      * Save the LinkMaster.
      *
      * \param ar output archive.
      */
      void saveLinkMaster(Serializable::OArchive& ar);
      #endif 

      #ifdef SIMP_TETHER
      /**
      * Read the TetherMaster.
      *
      * \param in input parameter stream
      */
      void readTetherMaster(std::istream& in);

      /**
      * Load the TetherMaster.
      *
      * \param ar input/loading archive
      */
      void loadTetherMaster(Serializable::IArchive& ar);

      /**
      * Save the TetherMaster.
      *
      * \param ar output/saving archive
      */
      void saveTetherMaster(Serializable::OArchive& ar);
      #endif // SIMP_TETHER

   private:

      /**
      * Pointer to DArray containing one MoleculeSet for each Species.
      *
      * MoleculeSet (*moleculeSetsPtr_)[i] contains all molecules in
      * this System that belong to Species i of the parent simulation.
      */
      DArray<MoleculeSet>* moleculeSetsPtr_;
 
      /// Pointer to Boundary object.
      Boundary* boundaryPtr_;
 
      #ifdef MCMD_LINK
      /// LinkMaster object to manage Links.
      LinkMaster* linkMasterPtr_;
      #endif

      #ifdef SIMP_TETHER
      /// TetherMaster object to manage Tethers.
      TetherMaster* tetherMasterPtr_;
      #endif

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      /// Pointer to the EnergyEnsemble.
      EnergyEnsemble* energyEnsemblePtr_;
   
      /// Pointer to the BoundaryEnsemble.
      BoundaryEnsemble* boundaryEnsemblePtr_;
  
      #ifndef SIMP_NOPAIR 
      /// Pointer to the PairPotential factory.
      PairFactory* pairFactoryPtr_;
      #endif
   
      #ifdef SIMP_BOND
      /// Pointer to the Factory<BondPotential>.
      Factory<BondPotential>* bondFactoryPtr_;
      #endif
  
      #ifdef SIMP_ANGLE 
      /// Pointer to the AnglePotential Factory.
      Factory<AnglePotential>* angleFactoryPtr_;
      #endif
   
      #ifdef SIMP_DIHEDRAL
      /// Pointer to DihedralPotential Factory
      Factory<DihedralPotential>* dihedralFactoryPtr_;
      #endif
  
      #ifdef SIMP_COULOMB
      /// Pointer to CoulombPotential Factory
      CoulombFactory*  coulombFactoryPtr_;
      #endif
  
      #ifdef SIMP_EXTERNAL
      /// Pointer to ExternalPotential factory
      Factory<ExternalPotential>* externalFactoryPtr_;
      #endif
  
      #ifdef MCMD_LINK
      /// Pointer to Link Factory
      Factory<BondPotential>* linkFactoryPtr_;
      #endif
   
      #ifdef SIMP_TETHER
      /// Pointer to TetherFactory.
      TetherFactory* tetherFactoryPtr_;
      #endif
   
      /// Pointer to a configuration reader/writer.
      ConfigIo* configIoPtr_;
   
      /// Pointer to a configuration reader/writer factory.
      Factory<ConfigIo>* configIoFactoryPtr_;
   
      /// Pointer to a trajectory reader/writer factory.
      Factory<TrajectoryReader>* trajectoryReaderFactoryPtr_;

      /// Pointer to a FileMaster.
      FileMaster* fileMasterPtr_;
   
      #ifdef MCMD_PERTURB
      /// Pointer to a perturbation object.
      Perturbation* perturbationPtr_;

      /// Pointer to a perturbation Factory.
      Factory<Perturbation>* perturbationFactoryPtr_;
      
      #ifdef UTIL_MPI
      /// Pointer to a ReplicaMove.
      ReplicaMove* replicaMovePtr_;
      
      /// Does system have ReplicaMove.
      bool hasReplicaMove_;

      #endif  // UTIL_MPI
      #endif  // MCMD_PERTURB
      #ifndef SIMP_NOPAIR
      /// Name of pair potential style.
      std::string pairStyle_;
      #endif

      #ifdef SIMP_BOND
      /// Name of bond potential style.
      std::string bondStyle_;
      #endif

      #ifdef SIMP_ANGLE
      /// Name of angle potential style.
      std::string angleStyle_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Name of dihedral potential style.
      std::string dihedralStyle_;
      #endif

      #ifdef SIMP_COULOMB
      /// Name of coulomb potential style.
      std::string coulombStyle_;
      #endif

      #ifdef SIMP_EXTERNAL
      /// Name of external potential style.
      std::string externalStyle_;
      #endif

      #ifdef MCMD_LINK
      /// Name of link potential style.
      std::string linkStyle_;
      #endif

      #ifdef SIMP_TETHER
      /// Name of tether potential style.
      std::string tetherStyle_;
      #endif

      /// Integer index for this System.
      int id_;

      /// Was this System instantiated with the copy constructor?
      bool isCopy_;

      /// Did this System instantiate a FileMaster object?
      bool createdFileMaster_;

      #ifdef MCMD_PERTURB
      /// Should this system read a Perturbation in the param file?
      bool expectPerturbationParam_;

      /// Has this System instantiated a Perturbation?
      bool createdPerturbation_;

      /// Has this System instantiated a PerturbationFactory?
      bool createdPerturbationFactory_;
      
      #ifdef UTIL_MPI
      /// Has this System instantiated a ReplicaMove?
      bool createdReplicaMove_;
      #endif // ifdef UTIL_MPI
      #endif // ifdef MCMD_PERTURB

      /// list of observers
      std::set<MoleculeSetObserver*> observers_;

      /// notify all observers
      void notifyMoleculeSetObservers() const;

   //friends:

      friend class SystemInterface;
      friend class ::SystemTest;

   }; 


   // Inline functions 

   /* 
   * Get integer Id for this System.
   */
   inline int System::id() const
   { return id_; }

   /* 
   * Get the parent Simulation by reference.
   */
   inline Simulation& System::simulation() const
   { 
      assert(simulationPtr_);
      return *simulationPtr_; 
   }
 
   /* 
   * Get the Boundary by reference.
   */
   inline Boundary& System::boundary() const
   { 
      assert(boundaryPtr_);
      return *boundaryPtr_;
   }

   #ifdef MCMD_LINK
   /* 
   * Get the LinkMaster by reference.
   */
   inline LinkMaster& System::linkMaster() const
   { 
      assert(linkMasterPtr_);
      return *linkMasterPtr_; 
   }
   #endif

   #ifdef SIMP_TETHER
   /* 
   * Get the TetherMaster by reference.
   */
   inline TetherMaster& System::tetherMaster() const
   { 
      assert(tetherMasterPtr_);
      return *tetherMasterPtr_; 
   }
   #endif

   /* 
   * Get the EnergyEnsemble by reference.
   */
   inline EnergyEnsemble& System::energyEnsemble() const
   { 
      assert(energyEnsemblePtr_);
      return *energyEnsemblePtr_; 
   }

   /* 
   * Get the BoundaryEnsemble by reference.
   */
   inline BoundaryEnsemble& System::boundaryEnsemble() const
   { 
      assert(boundaryEnsemblePtr_);
      return *boundaryEnsemblePtr_; 
   }

   /* 
   * Get the FileMaster by reference.
   */
   inline FileMaster& System::fileMaster() const
   { 
      assert(fileMasterPtr_);
      return *fileMasterPtr_; 
   }

   /* 
   * Was this System instantiated with the copy constructor?
   */
   inline bool System::isCopy() const
   { return isCopy_; }

   /* 
   * Get the number of molecules of a specific Species in this System.
   */
   inline int System::nMolecule(int speciesId) const
   {  
      assert(moleculeSetsPtr_);  
      return (*moleculeSetsPtr_)[speciesId].size(); 
   }

   /* 
   * Get a specific molecule of a specific Species.
   */
   inline Molecule& System::molecule(int speciesId, int moleculeId)
   {  
      assert(moleculeSetsPtr_);  
      return (*moleculeSetsPtr_)[speciesId][moleculeId]; 
   }

   /* 
   * Initialize a MoleculeIterator for molecules of one Species.
   */
   inline 
   void System::begin(int speciesId, MoleculeIterator& iterator)
   {
      assert(moleculeSetsPtr_);  
      (*moleculeSetsPtr_)[speciesId].begin(iterator); 
   }

   /* 
   * Initialize a MoleculeIterator for molecules of one Species.
   */
   inline 
   void System::begin(int speciesId, ConstMoleculeIterator& iterator) const
   {
      assert(moleculeSetsPtr_);  
      (*moleculeSetsPtr_)[speciesId].begin(iterator); 
   }

   #ifdef MCMD_PERTURB
   /**
   * Should this system create a Perturbation?
   */
   inline bool System::expectPerturbation() const
   { return expectPerturbationParam_; }

   /**
   * Return the perturbation factory by reference.
   */
   inline Factory<Perturbation>& System::perturbationFactory()
   { return *perturbationFactoryPtr_; }

   /* 
   * Does this system have an associated Perturbation?
   */
   inline bool System::hasPerturbation() const
   {  return perturbationPtr_; }

   /*
   * Get the Perturbation by reference.
   */
   inline Perturbation& System::perturbation() const
   { 
      assert(perturbationPtr_);
      return *perturbationPtr_; 
   }

   #ifdef UTIL_MPI
   /* 
   * Does this system have an associated ReplicaMove?
   */
   inline bool System::hasReplicaMove() const
   {  return replicaMovePtr_; }
   
   /*
   * Get the ReplicaMove by reference.
   */
   inline ReplicaMove& System::replicaMove() const
   {
      assert(replicaMovePtr_);
      return *replicaMovePtr_;
   }
   #endif // UTIL_MPI
   #endif // MCMD_PERTURB

   /*
   * Subscribe to moleculeSet change signal.
   */
   inline void System::subscribeMoleculeSetChange(MoleculeSetObserver& observer)
   {  observers_.insert(&observer); }

   /*
   * Unsubscribe from moleculeSet change signal.
   */
   inline void System::unsubscribeMoleculeSetChange(MoleculeSetObserver& observer)
   {  observers_.erase(&observer); }

   /*
   * Notifiy all observers
   */
   inline void System::notifyMoleculeSetObservers() const
   {
      std::set<MoleculeSetObserver*>::iterator itr;

      for ( itr = observers_.begin();
            itr != observers_.end(); itr++ )
         (*itr)->notifyMoleculeSetChanged();
   }

}
#endif
