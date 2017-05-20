#ifndef DDMD_SIMULATION_H
#define DDMD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class
#include <ddMd/communicate/Domain.h>             // member
#include <ddMd/communicate/Buffer.h>             // member
#include <ddMd/communicate/Exchanger.h>          // member
#include <ddMd/storage/AtomStorage.h>            // member
#include <ddMd/storage/BondStorage.h>            // member
#include <ddMd/storage/AngleStorage.h>           // member
#include <ddMd/storage/DihedralStorage.h>        // member
#include <ddMd/chemistry/AtomType.h>             // member (template param)
#include <ddMd/chemistry/MaskPolicy.h>           // member
#include <util/random/Random.h>                  // member
#include <util/boundary/Boundary.h>              // member
#include <util/space/Tensor.h>                   // member (template param)
#include <util/containers/DArray.h>              // member (template)
#include <util/containers/DMatrix.h>             // member (template)
#include <util/misc/Setable.h>                   // member (template)
#include <util/signal/Signal.h>                  // members
#include <util/archives/Serializable.h>          // typedef in function interface

#include <fstream>

namespace Util {
   class FileMaster;
   template <typename T> class Factory;
   class EnergyEnsemble;
   class BoundaryEnsemble;
}

namespace DdMd
{

   class Integrator;
   class ConfigIo;
   class SerializeConfigIo;
   #ifdef DDMD_MODIFIERS
   class ModifierManager;
   #endif
   class AnalyzerManager;
   class PairPotential;
   #ifdef SIMP_BOND
   class BondPotential;
   #endif
   #ifdef SIMP_ANGLE
   class AnglePotential;
   #endif
   #ifdef SIMP_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef SIMP_EXTERNAL
   class ExternalPotential;
   #endif

   using namespace Util;

   /**
   * Main object for a domain-decomposition MD simulation.
   *
   * A DdMd::Simulation contains and coordinates all the components of a
   * parallel MD simulation.
   *
   * Most member functions of this class should be called simultaneously
   * on all processors. This should be assumed to be the case unless the 
   * documentation says otherwise.  Member functions that should only be 
   * called on the master processor of the communicator associated with 
   * a spatial grid decomposition (i.e.,the communicator associated with 
   * the Domain member) are identified explicitly as such.
   *
   * \ingroup DdMd_Simulation_Module
   */
   class Simulation : public ParamComposite
   {

   public:

      using ParamComposite::readParam;
      using ParamComposite::load;

      // Lifetime

      #ifdef UTIL_MPI
      /**
      * Constructor.
      *
      * \param communicator MPI communicator for MD processors.
      */
      Simulation(MPI::Intracomm& communicator = MPI::COMM_WORLD);
      #else
      /**
      * Constructor.
      */
      Simulation();
      #endif

      /**
      * Destructor.
      */
      ~Simulation();

      /// \name Initialization
      //@{

      /**
      * Process command line options.
      *
      * Options:
      *
      *   -e  Enable echoing of the parameter file to the log file as it
      *       is read.
      *
      *   -s  nSystem [int]
      *       Enable multi system simulation, using different groups of
      *       processors for different systems. Integer argument nSystem
      *       is the number of systems. The rank of the communicator passed
      *       to the constructor must be an integer multiple of nSystem.
      *
      *   -p  filename [string]
      *       Specifies the name of a parameter file. If not specified here, 
      *       the parameter file for a single-system simulation may be read
      *       from standard input by default.
      *
      *   -r  filename [string]
      *       Restarts a simulation from a checkpoint file. The name of
      *       the checkpoint file is restartFilename + ".rst", i..e., it is
      *       obtained by adding a suffix ".rst" to the filename argument.
      *       Also sets the default command file name to filename + ".cmd".
      *
      *   -c  filename [string]
      *       Specifies the name of a command file. If not specified here,
      *       the command file name may be specified in the FileMaster block
      *       of the parameter file. 
      *
      * The -p (parameter file) and -r (restart) options are incompatible:
      * A simulation may either be initialized from a parameter file or a 
      * restart file, but not both. A command file is required in both cases.
      *
      * This function is passed two arguments that should contain the command 
      * line arguments that were passed to the main program:
      *
      * \param argc number of command line arguments
      * \param argv C-array of C-string argument strings
      */
      void setOptions(int argc, char* const * argv);

      /**
      * Read parameters from default parameter file.
      */
      virtual void readParam();

      /**
      * Read parameters.
      *
      * This version includes the opening and closing blocks.
      * Does nothing and returns if load() was called previously.
      */
      virtual void readParam(std::istream& in);

      /**
      * Read body of the parameter block, allocate memory and initialize.
      * 
      * This function is called by the readParam() methods, which simply
      * add code to read opening and closing brackets.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set flag to identify if reverse communication is enabled.
      *
      * Reverse communication is off by default. It is used in communication
      * patterns in which forces involving several atoms (nonbonded pairs or 
      * covalent groups) are calculated on one processor and the forces on 
      * particles that are not owned by that processor are then communicated 
      * to the processor that owns each particle. The default mode is for such 
      * forces to be calculated independently on each processor that owns one 
      * or more atoms in a group, to avoid the extra "reverse" communication
      * step. 
      *
      * \param reverseUpdateFlag true if reverse communication is enabled.
      */
      void setReverseUpdateFlag(bool reverseUpdateFlag);

      //@}
      /// \name Serialization (Restart files)
      //@{

      /**
      * Load internal state from a restart file.
      *
      * Call on all processors. This function opens an archive file with
      * a name given by concatentaing filename + ".rst" on the ioProcessor,
      * calls load(Serializable::OArchive& ), and closes the file. It also
      * sets the default command file name to filename + ".cmd".
      *
      * \param filename base filename (add suffixes ".rst" and ".cmd")
      */
      void load(const std::string& filename);

      /**
      * Load parameters from a restart archive.
      *
      * Call on all processors, but loads from archive only on ioProcessor.
      * This function loads both parameter information and the system
      * configuration.
      *
      * Do not call this directly. Instead call load(const std::string&)
      * or load(Serializable::IArchive& ).
      *
      * \param ar input archive (open for reading on ioProcessor).
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to a restart file.
      *
      * Call on all processors. Only writes from communicator ioProcessor.
      * This function opens an archive file with a name given by filename
      * + ".rst" on the ioProcessor, calls save(Serializable::OArchive& ),
      * and closes the file.
      */
      void save(const std::string& filename);

      /**
      * Save internal state to restart archive.
      *
      * Used to implement save(std::string).
      */
      virtual void save(Serializable::OArchive& ar);

      //@}
      /// \name Run Time Actions
      //@{

      /**
      * Read and execute commands from the default command file.
      *
      * This function calls readCommands(std::istream&) internally. The
      * current default command file is defined by the FileMaster. The 
      * file name for the default file may be specified as a command 
      * line option (option -c) or in the parameter file.
      */
      void readCommands();

      /**
      * Read and execute commands from a specific command file.
      *
      * \pre Simulation is initialized, by call to (read|load)Parameters.
      * \pre AtomStorage coordinates are scaled / generalized (not Cartesian)
      *
      * \param in command file
      */
      void readCommands(std::istream &in);

      /**
      * Set random velocities chosen from a Boltzmann distribution.
      *
      * \param temperature absolute temperature kT, in energy units.
      */
      void setBoltzmannVelocities(double temperature);

      /**
      * Subtract system center of mass velocity from atom velocities.
      *
      * Call on all processors. Does all-reduce internally. 
      *
      * \return drift velocity of system before subtraction.
      */
      Vector removeDriftVelocity();

      /**
      * Set forces for all atoms to zero.
      */
      void zeroForces();

      //@}
      /// \name Config File IO
      //@{

      /**
      * Read configuration file on master and distribute atoms.
      *
      * Upon return, all processors should have all atoms and groups,
      * and a full set of ghost atoms, but values for the atomic
      * forces are undefined.
      *
      * \pre AtomStorage is set for scaled / generalized coordinates
      * \post AtomStorage is set for scaled / generalized coordinates
      *
      * \param filename  name of input configuration file
      */
      void readConfig(const std::string& filename);

      /**
      * Write configuration file.
      *
      * \pre AtomStorage coordinates are generalized / scaled
      *
      * \param filename  name of output configuration file
      */
      void writeConfig(const std::string& filename);

      /**
      * Get the configuration file reader/writer factory by reference.
      */
      Factory<ConfigIo>& configIoFactory();

      /**
      * Create a new configuration file reader/writer.
      *
      * This method creates a new instance of the specified subclass
      * of ConfigIo, and retains a pointer to the new object.
      *
      * \param classname  name of desired ConfigIo subclass
      */
      void setConfigIo(std::string& classname);

      //@}
      /// \name Force, energy and stress calculators
      //@{

      /**
      * Compute forces for all local atoms.
      *
      * Upon return, forces are correct for all local atoms. Values
      * of the forces on ghost atoms are generally undefined.
      *
      * This method zeros all forces, adds forces from all potential
      * energies, and carries out reverse communication if required
      * (i.e., if reverseUpdateFlag is true).
      */
      void computeForces();

      /**
      * Compute forces for all local atoms and virial stress.
      *
      * Upon return, forces are correct for all local atoms and virial
      * stress values are set for all Potential objects. Values of the
      * forces on ghost atoms are undefined.
      *
      * This functions zeros all forces, adds forces from all potential
      * energies, and carries out reverse communication if required
      * (i.e., if reverseUpdateFlag is true). It also calculates and
      * sets each virial stress contribution.
      */
      void computeForcesAndVirial();

      /**
      * Compute total kinetic energy.
      *
      * This function is called on on all processors. The value is 
      * accessed by subsequently calling kineticEnergy() on the master.
      */
      void computeKineticEnergy();

      /**
      * Return precomputed total kinetic energy.
      *
      * Call only on master node.
      *
      * \return total kinetic energy for all nodes on master.
      */
      double kineticEnergy();

      /**
      * Mark kinetic energy as unknown.
      */
      void unsetKineticEnergy();

      /**
      * Calculate and store total potential energy on all processors.
      *
      * Reduce operation: Must be called on all nodes.
      */
      void computePotentialEnergies();

      /**
      * Return sum of precomputed total potential energies.
      *
      * Call only on master processor, after computePotentialEnergies.
      * Calls the energy() methods of each of the potential classes.
      *
      * \return total potential energy (only correct on master node).
      */
      double potentialEnergy() const;

      /**
      * Mark all potential energies as unknown.
      */
      void unsetPotentialEnergies();

      /**
      * Compute pair energies for each pair of atom types.
      *
      * Reduce operation: Must be called on all nodes.
      */
      void computePairEnergies();

      /**
      * Return precomputed pair energies.
      *
      * Call only on master processor, after computePairEnergies.
      *
      * \return total pair energies (only correct on master node).
      */
      DMatrix<double> pairEnergies() const;

      #ifdef SIMP_EXTERNAL
      /**
      * Compute pair energies for each pair of atom types.
      *
      * Reduce operation: Must be called on all nodes.
      */
      void computeExternalEnergy();

      /**
      * Return precomputed pair energies.
      *
      * Call only on master processor, after computePairEnergies.
      *
      * \return total pair energies (only correct on master node).
      */
      double externalEnergy() const;
      #endif

      /**
      * Calculate and store kinetic stress.
      *
      * Reduce operation: Must be called on all nodes.
      */
      void computeKineticStress();

      /**
      * Return total kinetic stress.
      *
      * Call only on master processor, after computeKineticStress.
      *
      * \return total kinetic stress (only correct on master node).
      */
      Tensor kineticStress() const;

      /**
      * Return total kinetic pressure.
      *
      * Call only on master processor, after computeKineticStress.
      *
      * \return total kinetic pressure only correct on master node).
      */
      double kineticPressure() const;

      /**
      * Mark kinetic stress as unknown.
      */
      void unsetKineticStress();

      /**
      * Calculate and store all virial stress contributions.
      *
      * Reduce operation: Must be called on all nodes. This calls
      * the computeStress() method of each potential class.
      */
      void computeVirialStress();

      /**
      * Return total virial stress.
      *
      * Call only on master processor, after computeVirialStress. This
      * method calls the stress() method of each potential class (pair,
      * bond, etc.) and adds and returns the sum.
      *
      * \return total virial stress (only correct on master node).
      */
      Tensor virialStress() const;

      /**
      * Return total virial pressure.
      *
      * Call only on master processor, after computeVirialStress.
      * Similar to virialStress, but returns average of diagonal
      * elements.
      *
      * \return total virial pressure (only correct on master node).
      */
      double virialPressure() const;

      /**
      * Mark all virial stress contributions as unknown.
      *
      * This method calls unsetStress for each potential energy.
      */
      void unsetVirialStress();

      //@}
      /// \name Potential Energy Classes (Objects, Style Strings and Factories)
      //@{

      /**
      * Get the PairPotential by const reference.
      */
      const PairPotential& pairPotential() const;

      /**
      * Get the PairPotential by reference.
      */
      PairPotential& pairPotential();

      /**
      * Return nonbonded pair style string.
      */
      std::string pairStyle() const;

      /**
      * Get the Factory<PairPotential> by reference.
      */
      Factory<PairPotential>& pairFactory();

      #ifdef SIMP_BOND
      /**
      * Get the PairPotential by const reference.
      */
      const BondPotential& bondPotential() const;

      /**
      * Get the PairPotential by reference.
      */
      BondPotential& bondPotential();

      /**
      * Return covalent bond style string.
      */
      std::string bondStyle() const;

      /**
      * Get the Factory<BondPotential> by reference.
      */
      Factory<BondPotential>& bondFactory();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get the AnglePotential by const reference.
      */
      const AnglePotential& anglePotential() const;

      /**
      * Get the AnglePotential by reference.
      */
      AnglePotential& anglePotential();

      /**
      * Return angle potential style string.
      */
      std::string angleStyle() const;

      /**
      * Get the AngleFactory by reference.
      */
      Factory<AnglePotential>& angleFactory();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get the DihedralPotential by const reference.
      */
      const DihedralPotential& dihedralPotential() const;

      /**
      * Get the DihedralPotential by reference.
      */
      DihedralPotential& dihedralPotential();

      /**
      * Return dihedral potential style string.
      */
      std::string dihedralStyle() const;

      /**
      * Get the associated Dihedral Factory by reference.
      */
      Factory<DihedralPotential>& dihedralFactory();
      #endif

      #ifdef SIMP_EXTERNAL
      /**
      * Get the ExternalPotential by reference.
      */
      const ExternalPotential& externalPotential() const;

      /**
      * Get the ExternalPotential by reference.
      */
      ExternalPotential& externalPotential();

      /**
      * Return external potential style string.
      */
      std::string externalStyle() const;

      /**
      * Get the associated External Factory by reference.
      */
      Factory<ExternalPotential>& externalFactory();
      #endif

      //@}
      /// \name Atom and Group Containers
      //@{

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& atomStorage();

      #ifdef SIMP_BOND
      /**
      * Get the BondStorage by reference.
      */
      BondStorage& bondStorage();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get the AngleStorage by reference.
      */
      AngleStorage& angleStorage();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get the DihedralStorage by reference.
      */
      DihedralStorage& dihedralStorage();
      #endif

      //@}
      /// \name Miscellaneous Accessors (return members by reference)
      //@{

      /**
      * Get the Domain by reference.
      */
      Domain& domain();

      /**
      * Get the Boundary by reference.
      */
      Boundary& boundary();

      /**
      * Get an AtomType descriptor by reference.
      *
      * \param i index of desired atom type.
      */
      AtomType& atomType(int i);

      /**
      * Get the Integrator by reference.
      */
      Integrator& integrator();

      /**
      * Get the EnergyEnsemble by reference.
      */
      EnergyEnsemble& energyEnsemble();

      /**
      * Get the BoundaryEnsemble by reference.
      */
      BoundaryEnsemble& boundaryEnsemble();

      /**
      * Get the associated FileMaster by reference.
      */
      FileMaster& fileMaster();

      /**
      * Get the Random number generator by reference.
      */
      Random& random();

      /**
      * Get the Exchanger by reference.
      */
      Exchanger& exchanger();

      /**
      * Get the Buffer by reference.
      */
      Buffer& buffer();

      #ifdef DDMD_MODIFIERS
      /**
      * Return the ModifierManager by reference.
      */
      ModifierManager& modifierManager();
      #endif

      /**
      * Return the AnalyzerManager by reference.
      */
      AnalyzerManager& analyzerManager();

      /**
      * Get the Integrator factory by reference.
      */
      Factory<Integrator>& integratorFactory();

      //@}
      /// \name Accessors (return by value)
      //@{

      /**
      * Get maximum number of atom types.
      */
      int nAtomType();

      #ifdef SIMP_BOND
      /**
      * Get maximum number of bond types.
      */
      int nBondType();
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get maximum number of angle types.
      */
      int nAngleType();
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get maximum number of dihedral types.
      */
      int nDihedralType();
      #endif

      #ifdef SIMP_EXTERNAL
      /**
      * Does this simulation have an external potential?
      */
      bool hasExternal();
      #endif

      /**
      * Return the value of the mask policy (MaskNone or MaskBonded).
      */
      MaskPolicy maskedPairPolicy() const;

      /**
      * Is reverse communication enabled?
      */
      bool reverseUpdateFlag() const;

      //@}
      /// \name Signals
      //@{

      /**
      * Signal to force unsetting of quantities that depend on x, v, or f.
      */
      Signal<>& modifySignal();

      /**
      * Signal to indicate change in atomic positions.
      */
      Signal<>& positionSignal();

      /**
      * Signal to indicate change in atomic velocities.
      */
      Signal<>& velocitySignal();

      /**
      * Signal to indicate exchange of atoms and groups.
      */
      Signal<>& exchangeSignal();

      //@}

      /**
      * Return true if this Simulation is valid, or throw an Exception.
      */
      bool isValid();

   protected:

      /**
      * Read the FileMaster.
      *
      * \param in input parameter stream
      */
      void readFileMaster(std::istream& in);

      /**
      * Load FileMaster from an archive.
      *
      * \param ar input/loading archive
      */
      void loadFileMaster(Serializable::IArchive& ar);

      /**
      * Save FileMaster to archive.
      *
      * \param ar output/saving archive
      */
      void saveFileMaster(Serializable::OArchive& ar);

      /**
      * Read potential styles and maskedPairPolicy.
      *
      * \param in input parameter stream
      */
      void readPotentialStyles(std::istream& in);

      /**
      * Load potential styles.
      *
      * \param ar input/loading archive
      */
      void loadPotentialStyles(Serializable::IArchive& ar);

      /**
      * Save potential styles to archive.
      *
      * \param ar output/saving archive
      */
      void savePotentialStyles(Serializable::OArchive& ar);

      /**
      * Read energy and boundary ensembles.
      *
      * \param in input parameter stream
      */
      void readEnsembles(std::istream& in);

      /**
      * Load energy and boundary ensembles from an input archive.
      *
      * \param ar input/loading archive
      */
      void loadEnsembles(Serializable::IArchive& ar);

      /**
      * Save energy and boundary ensembles to archive.
      *
      * \param ar output/saving archive
      */
      void saveEnsembles(Serializable::OArchive& ar);

   private:

      /// Container for all atoms and ghosts.
      AtomStorage atomStorage_;

      #ifdef SIMP_BOND
      /// Container for bonds.
      BondStorage bondStorage_;
      #endif

      #ifdef SIMP_ANGLE
      /// Container for angles.
      AngleStorage angleStorage_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Container for dihedrals.
      DihedralStorage dihedralStorage_;
      #endif

      /// Periodic system boundary.
      Boundary boundary_;

      /// Array of AtomType objects for all atoms in a simulation.
      DArray<AtomType> atomTypes_;

      /// Processor grid.
      Domain domain_;

      /// Communication buffer for sending atoms.
      Buffer buffer_;

      /// Exchanges atoms and ghosts for domain decomposition algorithm.
      Exchanger exchanger_;

      /// Random number generator.
      Random random_;

      /// Maximum boundary (used to allocate memory for the cell list).
      Boundary maxBoundary_;

      // Value of kinetic energy. Only valid on master.
      Setable<double> kineticEnergy_;

      // Value of kinetic stress. Only valid on master.
      Setable<Tensor> kineticStress_;

      /// Pointer to force/energy evaluator.
      PairPotential* pairPotentialPtr_;

      #ifdef SIMP_BOND
      /// Pointer to covalent bond potential.
      BondPotential* bondPotentialPtr_;
      #endif

      #ifdef SIMP_ANGLE
      /// Pointer to covalent 3-body angle potential.
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Pointer to covalent 4-body dihedral potential.
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef SIMP_EXTERNAL
      /// Pointer to external 1-body potential.
      ExternalPotential* externalPotentialPtr_;
      #endif

      /// Pointer to MD integrator.
      Integrator* integratorPtr_;

      /// Pointer to an EnergyEnsemble.
      EnergyEnsemble* energyEnsemblePtr_;

      /// Pointer to an BoundaryEnsemble.
      BoundaryEnsemble* boundaryEnsemblePtr_;

      /// Pointer to a FileMaster.
      FileMaster* fileMasterPtr_;

      /// Pointer to a configuration reader/writer.
      ConfigIo* configIoPtr_;

      /// Pointer to a configuration reader/writer for restart.
      SerializeConfigIo* serializeConfigIoPtr_;

      #ifdef DDMD_MODIFIERS
      /// ModifierManager
      ModifierManager* modifierManagerPtr_;
      #endif

      /// AnalyzerManager
      AnalyzerManager* analyzerManagerPtr_;

      /// Pointer to a PairPotential factory.
      Factory<PairPotential>* pairFactoryPtr_;

      #ifdef SIMP_BOND
      /// Pointer to a Factory<BondPotential>.
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

      #ifdef SIMP_EXTERNAL
      /// Pointer to ExternalPotential Factory
      Factory<ExternalPotential>* externalFactoryPtr_;
      #endif

      /// Pointer to MD integrator factory.
      Factory<Integrator>* integratorFactoryPtr_;

      /// Pointer to a configuration reader/writer factory.
      Factory<ConfigIo>* configIoFactoryPtr_;

      /// Name of pair potential style.
      std::string pairStyle_;

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

      #ifdef SIMP_EXTERNAL
      /// Name of external potential style.
      std::string externalStyle_;
      #endif

      /// Number of distinct atom types.
      int nAtomType_;

      #ifdef SIMP_BOND
      /// Number of distinct bond types.
      int nBondType_;
      #endif

      #ifdef SIMP_ANGLE
      /// Number of distinct angle types.
      int nAngleType_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Number of distinct dihedral types.
      int nDihedralType_;
      #endif

      #ifdef SIMP_EXTERNAL
      /// Does this simulation have an external potential?
      bool hasExternal_;
      #endif

      /// Does this simulation keep track of AtomContext info?
      bool hasAtomContext_;

      /**
      * Policy for suppressing pair interactions for some atom pairs.
      *
      * Allowed values of enum MaskPolicy:
      *
      *  - MaskNone:   no masked pairs
      *  - MaskBonded:  mask pair interaction between bonded atoms
      */
      MaskPolicy maskedPairPolicy_;

      /// Is reverse communication enabled?
      bool reverseUpdateFlag_;

      #ifdef UTIL_MPI
      /// Communicator for this system.
      MPI::Intracomm communicator_;
      #endif

      /// Signal to force clearing of all computed quantities.
      Signal<>  modifySignal_;

      /// Signal to indicate change in atomic positions.
      Signal<>  positionSignal_;

      /// Signal to indicate change in atomic velocities.
      Signal<>  velocitySignal_;

      /// Signal to indicate exchange of atoms ownership.
      Signal<>  exchangeSignal_;

      /// Log output file (if not standard out)
      std::ofstream logFile_;

      /// Has readParam or load been called?
      bool isInitialized_;

      /// Is this Simulation in the process of restarting?
      bool isRestarting_;

      /// Return the current ConfigIo (create if necessary)
      ConfigIo& configIo();

      /// Return a SerializeConfigIo (create if necessary)
      SerializeConfigIo& serializeConfigIo();

      void setGroup(std::stringstream& inBuffer);

   // friends:

      friend class SimulationAccess;

   };

   // Inline method definitions

   inline Boundary& Simulation::boundary()
   { return boundary_; }

   inline Domain& Simulation::domain()
   { return domain_; }

   inline AtomStorage& Simulation::atomStorage()
   { return atomStorage_; }

   #ifdef SIMP_BOND
   inline BondStorage& Simulation::bondStorage()
   { return bondStorage_; }
   #endif

   #ifdef SIMP_ANGLE
   inline AngleStorage& Simulation::angleStorage()
   { return angleStorage_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   inline DihedralStorage& Simulation::dihedralStorage()
   { return dihedralStorage_; }
   #endif

   inline Exchanger& Simulation::exchanger()
   { return exchanger_; }

   inline Buffer& Simulation::buffer()
   { return buffer_; }

   inline const PairPotential& Simulation::pairPotential() const
   {
      assert(pairPotentialPtr_);
      return *pairPotentialPtr_;
   }

   inline PairPotential& Simulation::pairPotential()
   {
      assert(pairPotentialPtr_);
      return *pairPotentialPtr_;
   }

   #ifdef SIMP_BOND
   inline const BondPotential& Simulation::bondPotential() const
   {
      assert(bondPotentialPtr_);
      return *bondPotentialPtr_;
   }

   inline BondPotential& Simulation::bondPotential()
   {
      assert(bondPotentialPtr_);
      return *bondPotentialPtr_;
   }
   #endif

   #ifdef SIMP_ANGLE
   inline const AnglePotential& Simulation::anglePotential() const
   {
      assert(anglePotentialPtr_);
      return *anglePotentialPtr_;
   }

   inline AnglePotential& Simulation::anglePotential()
   {
      assert(anglePotentialPtr_);
      return *anglePotentialPtr_;
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   inline const DihedralPotential& Simulation::dihedralPotential() const
   {
      assert(dihedralPotentialPtr_);
      return *dihedralPotentialPtr_;
   }

   inline DihedralPotential& Simulation::dihedralPotential()
   {
      assert(dihedralPotentialPtr_);
      return *dihedralPotentialPtr_;
   }
   #endif

   #ifdef SIMP_EXTERNAL
   inline const ExternalPotential& Simulation::externalPotential() const
   {
      assert(externalPotentialPtr_);
      return *externalPotentialPtr_;
   }

   inline ExternalPotential& Simulation::externalPotential()
   {
      assert(externalPotentialPtr_);
      return *externalPotentialPtr_;
   }
   #endif

   /// Get the MD integrator by reference.
   inline Integrator& Simulation::integrator()
   {
      assert(integratorPtr_);
      return *integratorPtr_;
   }

   /// Get the Random RNG by reference.
   inline Random& Simulation::random()
   { return random_; }

   /// Get the EnergyEnsemble by reference.
   inline EnergyEnsemble& Simulation::energyEnsemble()
   {
      assert(energyEnsemblePtr_);
      return *energyEnsemblePtr_;
   }

   /// Get the BoundaryEnsemble by reference.
   inline BoundaryEnsemble& Simulation::boundaryEnsemble()
   {
      assert(boundaryEnsemblePtr_);
      return *boundaryEnsemblePtr_;
   }

   /// Get the FileMaster by reference.
   inline FileMaster& Simulation::fileMaster()
   {
      assert(fileMasterPtr_);
      return *fileMasterPtr_;
   }

   #ifdef DDMD_MODIFIERS
   /// Get the ModifierManager by reference.
   inline ModifierManager& Simulation::modifierManager()
   {
      assert(modifierManagerPtr_);
      return *modifierManagerPtr_;
   }
   #endif

   /// Get the AnalyzerManager by reference.
   inline AnalyzerManager& Simulation::analyzerManager()
   {
      assert(analyzerManagerPtr_);
      return *analyzerManagerPtr_;
   }

   /// Get maximum number of atom types.
   inline int Simulation::nAtomType()
   {  return nAtomType_; }

   #ifdef SIMP_BOND
   /// Get maximum number of bond types.
   inline int Simulation::nBondType()
   {  return nBondType_; }
   #endif

   #ifdef SIMP_ANGLE
   /// Get maximum number of angle types.
   inline int Simulation::nAngleType()
   {  return nAngleType_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /// Get the maximum number of dihedral types.
   inline int Simulation::nDihedralType()
   {  return nDihedralType_; }
   #endif

   #ifdef SIMP_EXTERNAL
   /// Does this simulation have an external potential?
   inline bool Simulation::hasExternal()
   {  return hasExternal_; }
   #endif

   /// Get an AtomType descriptor for a specific type by reference.
   inline AtomType& Simulation::atomType(int i)
   {  return atomTypes_[i]; }

   /// Get the masked pair policy.
   inline MaskPolicy Simulation::maskedPairPolicy() const
   {  return maskedPairPolicy_; }

   /// Get the reverse update flag.
   inline bool Simulation::reverseUpdateFlag() const
   {  return reverseUpdateFlag_; }

   // Signals

   /// Signal to force unsetting of all computed quantities.
   inline Signal<>& Simulation::modifySignal()
   {  return modifySignal_; }

   /// Signal to indicate change in atomic positions.
   inline Signal<>& Simulation::positionSignal()
   { return positionSignal_; }

   /// Signal to indicate change in atomic velocities.
   inline Signal<>& Simulation::velocitySignal()
   { return velocitySignal_; }

   /// Signal to indicate exchange of atom ownership.
   inline Signal<>& Simulation::exchangeSignal()
   { return exchangeSignal_; }

}
#endif
