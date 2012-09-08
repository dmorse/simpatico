#ifndef MCMD_SYSTEM_H
#define MCMD_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <util/param/ParamComposite.h>        // base class
#include <util/archives/Serializable.h>       // base class

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
   class FileMaster;
   class ConfigIo;
   class TrajectoryIo;
   class PairFactory;
   class BondPotential;
   #ifdef INTER_ANGLE
   class AnglePotential;
   #endif
   #ifdef INTER_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef MCMD_LINK
   class LinkPotential;
   class LinkMaster;
   #endif
   #ifdef INTER_EXTERNAL
   class ExternalPotential;
   #endif
   #ifdef INTER_TETHER
   class TetherFactory;
   class TetherMaster;
   #endif
   #ifdef MCMD_PERTURB
   class Perturbation;
   class ReplicaMove;
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
   class System : public ParamComposite, public Serializable
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
      * Can be used to set the FileMaster to that of the parent 
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
      * copy, this method does nothing and returns normally.
      *
      * \param in pararameter file input stream
      */
      virtual void readParameters(std::istream& in);
 
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
      * This method allows one to choose from among several subclasses
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
      * This method uses a ConfigIo object that is registered
      * with this System. If a ConfigIo object has not been
      * registered by calling setConfigIo(std::string&), this 
      * method creates and uses an instance of the ConfigIo 
      * base class.
      *
      * Precondition: The System and its Parent Simulation must
      * have been initialized by calling their readParam methods. 
      *
      * \param in configuration file input stream
      */
      virtual void readConfig(std::istream& in);

      /**
      * Write system configuration to a specified ostream.
      *
      * Like readConfig(), this method will create and use 
      * a ConfigIo object if none has been registered 
      * previously. 
      *
      * \param out configuration file output stream
      */
      void writeConfig(std::ostream& out);

      /**
      * Serialize the System to/from an archive.
      *
      * \param ar      output or input Archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Save the System configuration to an archive.
      *
      * \param ar output (saving) archive object.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load the System configuration from an archive.
      *
      * \param ar input (loading) archive object.
      */
      virtual void load(Serializable::IArchiveType& ar);

      //@}
      /// \name Trajectory File IO
      //@{

      /**
      * Get the trajectory reader/writer factory by reference.
      */
      Factory<TrajectoryIo>& trajectoryIoFactory();

      //@}
      /// \name Molecule Set Mutators
      //@{
      
      /**
      * Add a Molecule to this System.  
      *
      * This method adds a Molecule to the set of Molecules of the same
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
      * This method removes a Molecule from the set of molecules of the
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
      * This method returns the current (mutable) index of the molecule 
      * within the set of molecules of the same Species in this System, 
      * in the range 0 <= moleculeId < nMolecule(speciesId). This is 
      * the index required as the second argument of molecule(int, int).
      *
      * \param molecule Molecule object of interest.
      * \return index for molecule within its Species and System.
      */
      int moleculeId(const Molecule& molecule) const;

      /** 
      * Get a specific Molecule in this System, by integer index.
      *
      * The moleculeId must be in range 0 <= moleculeId < nMolecule(speciesId).
      * The index associated with a molecule is volatile, and can change when
      * another molecule of the same Species is removed from this System. 
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
    
      #ifndef INTER_NOPAIR 

      /**
      * Get the PairFactory by reference.
      */
      PairFactory& pairFactory();

      /**
      * Return nonbonded pair style string.
      */
      std::string pairStyle() const;
  
      #endif

      /**
      * Get the associated Factory<BondPotential> by reference.
      */
      Factory<BondPotential>& bondFactory();

      /**
      * Return covalent bond style string.
      */
      std::string bondStyle() const;

      #ifdef INTER_ANGLE
      /**
      * Get the associated AngleFactory by reference.
      */
      Factory<AnglePotential>& angleFactory();

      /**
      * Return angle potential style string.
      */
      std::string angleStyle() const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get the associated Dihedral Factory by reference.
      */
      Factory<DihedralPotential>& dihedralFactory();

      /**
      * Return dihedral potential style string.
      */
      std::string dihedralStyle() const;
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

      #ifdef INTER_EXTERNAL
      /**
      * Get the associated ExternalPotential factory by reference.
      */
      Factory<ExternalPotential>& externalFactory();

      /**
      * Return external potential style string.
      */
      std::string externalStyle() const;
      #endif

      #ifdef INTER_TETHER
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

      //@}
      #endif

      #endif

      /// \name Accessors (Miscellaneous)
      //@{
     
      /// Get integer Id for this System.
      int id() const;
  
      /// Get the parent Simulation by reference.
      Simulation& simulation() const;

      /// Get the Boundary by reference.
      Boundary& boundary() const;

      /// Get the EnergyEnsemble by reference.
      EnergyEnsemble& energyEnsemble() const;

      /// Get the BoundaryEnsemble by reference.
      BoundaryEnsemble& boundaryEnsemble() const;

      /// Get the associated FileMaster by reference.
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
      * Return a pointer to a new default TrajectoryIoFactory.
      */
      virtual Factory<TrajectoryIo>* newDefaultTrajectoryIoFactory();

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

      #ifdef UTIL_MPI
      /**
      * Read the ReplicaMove parameter block (if any)
      *
      * \param in input parameter stream
      */
      void readReplicaMove(std::istream& in);
      #endif

      #endif
      /**
      * Allocate and initialize molecule sets for all species.
      *
      * This method is called within Simulation::initialize() private 
      * method to allocate and initialize an array of MoleculeSet objects 
      * for all Species for this System.
      *
      * Preconditions: This System must be associated with a parent
      * Simulation, and all of the Species objects must have been
      * initialized by calling SpeciesManager::readParameters().
      */
      void allocateMoleculeSets();

      /**
      * If no FileMaster exists, create and initialize one. 
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readFileMaster(std::istream& in);

      /**
      * Read potential styles, initialize LinkMaster or TetherMaster if needed.
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readPotentialStyles(std::istream& in);

      /**
      * Read energy and boundary ensembles.
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readEnsembles(std::istream& in);

      #ifdef MCMD_LINK
      /**
      * Read the LinkMaster.
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readLinkMaster(std::istream& in);
      #endif

      #ifdef INTER_TETHER
      /**
      * Read the TetherMaster.
      *
      * Invoked in implementation of readParameters().
      *
      * \param in input parameter stream
      */
      void readTetherMaster(std::istream& in);
      #endif

   private:

      /**
      * Pointer to DArray containing one MoleculeSet for each Species.
      *
      * The MoleculeSet (*moleculeSetsPtr_)[i] contains the molecules in
      * this System that belong to Species i of the parent simulation.
      */
      DArray< MoleculeSet >* moleculeSetsPtr_;
  
      /// Pointer to Boundary object for actual boundary.
      Boundary*         boundaryPtr_;
 
      #ifdef MCMD_LINK
      /// TetherMaster object to manage Tethers
      LinkMaster*       linkMasterPtr_;
      #endif

      #ifdef INTER_TETHER
      /// TetherMaster object to manage Tethers
      TetherMaster*     tetherMasterPtr_;
      #endif

      /// Pointer to parent Simulation.
      Simulation*       simulationPtr_;

      /// Pointer to an EnergyEnsemble.
      EnergyEnsemble*   energyEnsemblePtr_;
   
      /// Pointer to an BoundaryEnsemble.
      BoundaryEnsemble* boundaryEnsemblePtr_;
  
      #ifndef INTER_NOPAIR 
      /// Pointer to a PairPotential factory.
      PairFactory*  pairFactoryPtr_;
      #endif
   
      /// Pointer to a Factory<BondPotential>.
      Factory<BondPotential>*  bondFactoryPtr_;
  
      #ifdef INTER_ANGLE 
      /// Pointer to the AnglePotential Factory.
      Factory<AnglePotential>*  angleFactoryPtr_;
      #endif
   
      #ifdef INTER_DIHEDRAL
      /// Pointer to DihedralPotential Factory
      Factory<DihedralPotential>*  dihedralFactoryPtr_;
      #endif
  
      #ifdef MCMD_LINK
      /// Pointer to Link Factory
      Factory<BondPotential>*  linkFactoryPtr_;
      #endif
   
      #ifdef INTER_EXTERNAL
      /// Pointer to ExternalPotential factory
      Factory<ExternalPotential>*  externalFactoryPtr_;
      #endif
  
      #ifdef INTER_TETHER
      /// Pointer to TetherFactory
      TetherFactory*  tetherFactoryPtr_;
      #endif
   
      /// Pointer to a configuration reader/writer.
      ConfigIo*         configIoPtr_;
   
      /// Pointer to a configuration reader/writer factory.
      Factory<ConfigIo>* configIoFactoryPtr_;
   
      /// Pointer to a trajectory reader/writer factory.
      Factory<TrajectoryIo>* trajectoryIoFactoryPtr_;

      /// Pointer to a FileMaster.
      FileMaster*       fileMasterPtr_;
   
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

      #endif
      #endif 

      #ifndef INTER_NOPAIR
      /// Name of pair potential style.
      std::string pairStyle_;
      #endif

      /// Name of bond potential style.
      std::string bondStyle_;

      #ifdef INTER_ANGLE
      /// Name of angle potential style.
      std::string angleStyle_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Name of dihedral potential style.
      std::string dihedralStyle_;
      #endif

      #ifdef MCMD_LINK
      /// Name of link potential style.
      std::string linkStyle_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Name of external potential style.
      std::string externalStyle_;
      #endif

      #ifdef INTER_TETHER
      /// Name of tether potential style.
      std::string tetherStyle_;
      #endif

      /// Integer index for this System.
      int     id_;

      /// Was this System instantiated with the copy constructor?
      bool    isCopy_;

      /// Did this System instantiate a FileMaster object?
      bool    createdFileMaster_;

      #ifdef MCMD_PERTURB
      /// Should this system read a Perturbation in the param file?
      bool    expectPerturbationParam_;

      /// Has this System instantiated a Perturbation?
      bool    createdPerturbation_;

      /// Has this System instantiated a PerturbationFactory?
      bool    createdPerturbationFactory_;
      
      /// Has this System instantiated a ReplicaMove?
      bool    createdReplicaMove_;
      #endif


      /// list of observers
      std::set<MoleculeSetObserver*> observers_;

      /// notify all observers
      void    notifyMoleculeSetObservers() const;

   //friends:


      friend class SubSystem;

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

   #ifdef INTER_TETHER
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
   * return the perturbation factory by reference.
   */
   inline Factory<Perturbation>& System::perturbationFactory()
   { return *perturbationFactoryPtr_; }


   /* 
   * Does this system have an associated Perturbation?
   */
   inline bool System::hasPerturbation() const
   {  return (perturbationPtr_ != 0); }

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
   {  return (replicaMovePtr_ != 0); }
   
   /*
   * Get the ReplicaMove by reference.
   */
   inline ReplicaMove& System::replicaMove() const
   {
      assert(replicaMovePtr_);
      return *replicaMovePtr_;
   }
   #endif

   inline bool System::expectPerturbation() const
   { return expectPerturbationParam_; }

   #endif

   /*
    * Subscribe to moleculeSet change signal
    */
   inline void System::subscribeMoleculeSetChange(MoleculeSetObserver& observer)
   {
      observers_.insert(&observer);
   }

   /*
    * Unsubscribe from moleculeSet change signal
    */
   inline void System::unsubscribeMoleculeSetChange(MoleculeSetObserver& observer)
   {
      observers_.erase(&observer);
   }

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
