#ifndef DDMD_SIMULATION_H
#define DDMD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
#include <ddMd/chemistry/AtomType.h>             // member
#include <ddMd/chemistry/MaskPolicy.h>           // member
#include <ddMd/diagnostics/DiagnosticManager.h>  // member
#include <util/boundary/Boundary.h>              // member 
#include <util/random/Random.h>                  // member 
#include <util/containers/DArray.h>              // member 

namespace Util { 
   template <typename T> class Factory; 
   class EnergyEnsemble;
   class BoundaryEnsemble;
}
namespace McMd { 
   class McSimulation; 
}

namespace DdMd
{

   class ConfigIo;
   class FileMaster;
   class DiagnosticManager;
   class Integrator;
   class PairPotential;
   class BondPotential;
   #ifdef INTER_ANGLE
   class AnglePotential;
   #endif
   #ifdef INTER_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef INTER_EXTERNAL
   class ExternalPotential;
   #endif

   using namespace Util;

   /**
   * Main object for a domain-decomposition MD simulation.
   *
   * A DdMd::Simulation contains and coordinates all the components of a 
   * parallel MD simulation. 
   *
   * \ingroup DdMd_Simulation_Module
   */
   class Simulation : public ParamComposite
   {

   public:

      #ifdef UTIL_MPI

      /**
      * Default constructor.
      *
      * \param communicator MPI communicator for MD processors.
      */
      Simulation(MPI::Intracomm& communicator = MPI::COMM_WORLD);

      /**
      * Constructor.
      *
      * \param mcSimulation parent McSimulation object.
      * \param communicator MPI communicator for MD processors.
      */
      Simulation(McMd::McSimulation& mcSimulation,
             MPI::Intracomm& communicator = MPI::COMM_WORLD);
      #else

      /**
      * Default constructor.
      */
      Simulation();

      #endif

      /**
      * Destructor.
      */
      ~Simulation();

      /**
      * Process command line options.
      *
      * Options:
      *
      *   -e  Enable echoing of the parameter file to the log file as it
      *       is read.
      *
      * \param argc number of arguments
      * \param argv vector of C-string argument strings
      */
      void setOptions(int argc, char **argv);

      /**
      * Read parameters, allocate memory and initialize.
      */
      virtual void readParam(std::istream& in);

      /*
      * Read and execute commands from a command file.
      */
      void readCommands(std::istream &in);
   
      /**
      * Read and execute commands from the default command file.
      */
      void readCommands();

      // Mutators

      /**
      * Integrate equations of motion. 
      */
      void simulate(int nStep);

      /**
      * Set random velocities chosen from Boltzmann distribution.
      *  
      * \param temperature absolute temperature kT, in energy units. 
      */
      void setBoltzmannVelocities(double temperature);

      /**
      * Set forces for all local atoms to zero.
      */
      void zeroForces();

      /**
      * Set flag to identify if force communication is enabled.
      *
      * \param forceCommFlag true if force communication is enabled.
      */
      void setForceCommFlag(bool forceCommFlag);

      /// \name Config File IO
      //@{

      /**
      * Read configuration file on master and distribute atoms.
      *
      * Upon return, all processors should have all atoms and
      * groups, and a full set of ghost atoms, but values for
      * the atomic forces are undefined.
      *
      * \param filename name of configuration file.
      */
      void readConfig(const std::string& filename);

      /**
      * Write configuration file.
      *
      * \param filename name of configuration file.
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
      * \param classname name of desired ConfigIo subclass.
      */
      void setConfigIo(std::string& classname);

      //@}
      /// \name Force and Energy calculators
      //@{
      
      /**
      * Compute forces for all local atoms.
      *
      * Upon return, forces are correct for all local atoms. Values of
      * the forces on ghost atoms are undefined.
      *
      * This method zeros all forces, adds forces from all potential
      * energies, and carries out reverse force communication only if
      * required (i.e., if forceCommFlag is true). 
      */
      void computeForces();

      /**
      * Compute total kinetic energy.
      * 
      * Reduce operation: Must be called on all processors.
      * Value is returned by calling kineticEnergy on master.
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
      * Calculate and store total potential energy on all processors.
      * 
      * Reduce operation: Must be called on all nodes.
      *
      * \return total pair potential for all nodes on master, 0.0 otherwise.
      */
      void computePotentialEnergies();

      /**
      * Return precomputed total potential energy.
      *
      * Call only on master processor, after computePotentialEnergies.
      * 
      * \return total potential energy (only correct on master node).
      */
      double potentialEnergy();

      //@}
      /// \name Potential Energy Factories and Styles
      //@{

      #ifndef INTER_NOPAIR
      /**
      * Get the Factory<PairPotential> by reference.
      */
      Factory<PairPotential>& pairFactory();

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

      #ifdef INTER_EXTERNAL
      /**
      * Get the associated External Factory by reference.
      */
      Factory<ExternalPotential>& externalFactory();

      /**
      * Return external potential style string.
      */
      std::string externalStyle() const;
      #endif

      //@}
      /// \name Accessors (return members by reference)
      //@{

      /**
      * Get the Domain by reference.
      */
      Domain& domain();

      /**
      * Get the AtomStorage by reference.
      */
      AtomStorage& atomStorage();
   
      /**
      * Get the BondStorage by reference.
      */
      BondStorage& bondStorage();
  
      #ifdef INTER_ANGLE 
      /**
      * Get the angleStorage by reference.
      */
      AngleStorage& angleStorage();
      #endif
   
      #ifdef INTER_DIHEDRAL
      /**
      * Get the angleStorage by reference.
      */
      DihedralStorage& dihedralStorage();
      #endif
   
      /**
      * Get the Boundary by reference.
      */
      Boundary& boundary();
   
      /**
      * Get the PairPotential by reference.
      */
      PairPotential& pairPotential();
   
      /**
      * Get the PairPotential by reference.
      */
      BondPotential& bondPotential();
  
      #ifdef INTER_ANGLE 
      /**
      * Get the AnglePotential by reference.
      */
      AnglePotential& anglePotential();
      #endif
   
      #ifdef INTER_DIHEDRAL
      /**
      * Get the DihedralPotential by reference.
      */
      DihedralPotential& dihedralPotential();
      #endif
   
      #ifdef INTER_EXTERNAL
      /**
      * Get the ExternalPotential by reference.
      */
      ExternalPotential& externalPotential();
      #endif
   
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
   
      /// Get the EnergyEnsemble by reference.
      EnergyEnsemble& energyEnsemble();

      /// Get the BoundaryEnsemble by reference.
      BoundaryEnsemble& boundaryEnsemble();

      /// Get the associated FileMaster by reference.
      FileMaster& fileMaster();

      /**
      * Get the Md  integrator factory by reference.
      */
      Factory<Integrator>& integratorFactory();

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
  
      /**
      * Return associated DiagnosticManager by reference.
      */
      DiagnosticManager& diagnosticManager();

      //@}
      /// \name Accessors (return by value)
      //@{
      
      /**
      * Get maximum number of atom types.
      */
      int nAtomType();

      /**
      * Get maximum number of bond types.
      */
      int nBondType();

      #ifdef INTER_ANGLE
      /**
      * Get maximum number of angle types.
      */
      int nAngleType();
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get maximum number of dihedral types.
      */
      int nDihedralType();
      #endif

      #ifdef INTER_EXTERNAL
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
      * Is reverse force communication enabled?
      */
      bool forceCommFlag() const;

      //@}

      /**
      * Return true if this Simulation is valid, or throw an Exception.
      */
      bool isValid();

   protected:

      #if 0
      /**
      * Return a pointer to a new default ConfigIo.
      */
      virtual ConfigIo* newDefaultConfigIo();

      /**
      * Return a pointer to a new default ConfigIoFactory.
      */
      virtual Factory<ConfigIo>* newDefaultConfigIoFactory();
      #endif

      /**
      * Read the FileMaster.
      *
      * \param in input parameter stream
      */
      void readFileMaster(std::istream& in);

      /**
      * Read potential styles and maskedPairPolicy.
      *
      * \param in input parameter stream
      */
      void readPotentialStyles(std::istream& in);

      /**
      * Read energy and boundary ensembles.
      *
      * \param in input parameter stream
      */
      void readEnsembles(std::istream& in);

   private:

      /// Container for all atoms and ghosts.
      AtomStorage   atomStorage_;

      /// Container for bonds.
      BondStorage   bondStorage_;

      #ifdef INTER_ANGLE
      /// Container for angles.
      AngleStorage  angleStorage_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Container for dihedrals.
      DihedralStorage  dihedralStorage_;
      #endif

      /// Periodic system boundary.
      Boundary      boundary_;

      /// Array of AtomType objects for all atoms in a simulation.
      DArray<AtomType> atomTypes_;

      /// Processor grid.
      Domain        domain_;

      /// Communication buffer for sending atoms.
      Buffer        buffer_;

      /// Exchanges atoms and ghosts for domain decomposition algorithm.
      Exchanger     exchanger_;

      /// Random number generator.
      Random        random_;

      /// Maximum boundary (used to allocate memory for the cell list).
      Boundary      maxBoundary_;

      // Value of total kinetic energy, for all processors.
      double  kineticEnergy_;

      /// Pointer to force/energy evaluator.
      PairPotential* pairPotentialPtr_;

      /// Pointer to covalent bond potential.
      BondPotential* bondPotentialPtr_;

      #ifdef INTER_ANGLE
      /// Pointer to covalent 3-body angle potential.
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Pointer to covalent 4-body dihedral potential.
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Pointer to external 1-body potential.
      ExternalPotential* externalPotentialPtr_;
      #endif

      /// Pointer to MD integrator.
      Integrator*  integratorPtr_;

      /// Pointer to an EnergyEnsemble.
      EnergyEnsemble*  energyEnsemblePtr_;

      /// Pointer to an BoundaryEnsemble.
      BoundaryEnsemble*  boundaryEnsemblePtr_;

      /// Pointer to a FileMaster.
      FileMaster*         fileMasterPtr_;

      /// Pointer to a configuration reader/writer.
      ConfigIo*         configIoPtr_;
   
      /// DiagnosticManager
      DiagnosticManager*  diagnosticManagerPtr_;

      #ifndef INTER_NOPAIR
      /// Pointer to a PairPotential factory.
      Factory<PairPotential>*  pairFactoryPtr_;
      #endif

      /// Pointer to a Factory<BondPotential>.
      Factory<BondPotential>* bondFactoryPtr_;

      #ifdef INTER_ANGLE
      /// Pointer to the AnglePotential Factory.
      Factory<AnglePotential>* angleFactoryPtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Pointer to DihedralPotential Factory
      Factory<DihedralPotential>* dihedralFactoryPtr_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Pointer to ExternalPotential Factory
      Factory<ExternalPotential>* externalFactoryPtr_;
      #endif

      /// Pointer to MD integrator factory.
      Factory<Integrator>*  integratorFactoryPtr_;

      /// Pointer to a configuration reader/writer factory.
      Factory<ConfigIo>* configIoFactoryPtr_;

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

      #ifdef INTER_EXTERNAL
      /// Name of external potential style.
      std::string externalStyle_;
      #endif

      /// Number of distinct atom types.
      int nAtomType_;

      /// Number of distinct bond types.
      int nBondType_;

      #ifdef INTER_ANGLE
      /// Number of distinct angle types.
      int nAngleType_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Number of distinct dihedral types.
      int nDihedralType_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Does this simulation have an external potential?
      bool hasExternal_;
      #endif

      /**
      * Policy for suppressing pair interactions for some atom pairs.
      *
      * Allowed values of enum MaskPolicy:
      *
      *  - MaskNone:   no masked pairs
      *  - MaskBonded:  mask pair interaction between bonded atoms
      */
      MaskPolicy  maskedPairPolicy_;

      /// Is reverse force communication enabled?
      bool forceCommFlag_;

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

   inline BondStorage& Simulation::bondStorage()
   { return bondStorage_; }

   #ifdef INTER_ANGLE
   inline AngleStorage& Simulation::angleStorage()
   { return angleStorage_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline DihedralStorage& Simulation::dihedralStorage()
   { return dihedralStorage_; }
   #endif

   inline Exchanger& Simulation::exchanger()
   { return exchanger_; }

   inline Buffer& Simulation::buffer()
   { return buffer_; }

   inline PairPotential& Simulation::pairPotential()
   { 
      assert(pairPotentialPtr_); 
      return *pairPotentialPtr_; 
   }

   inline BondPotential& Simulation::bondPotential()
   { 
      assert(bondPotentialPtr_); 
      return *bondPotentialPtr_; 
   }

   #ifdef INTER_ANGLE
   inline AnglePotential& Simulation::anglePotential()
   { 
      assert(anglePotentialPtr_); 
      return *anglePotentialPtr_; 
   }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Get the DihedralPotential by reference.
   */
   inline DihedralPotential& Simulation::dihedralPotential()
   { 
      assert(dihedralPotentialPtr_); 
      return *dihedralPotentialPtr_; 
   }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Get the ExternalPotential by reference.
   */
   inline ExternalPotential& Simulation::externalPotential()
   { 
      assert(externalPotentialPtr_); 
      return *externalPotentialPtr_; 
   }
   #endif

   /*
   * Get the Integrator by reference.
   */
   inline Integrator& Simulation::integrator()
   {
      assert(integratorPtr_); 
      return *integratorPtr_; 
   }

   inline Random& Simulation::random()
   { return random_; }

   /*
   * Get the EnergyEnsemble by reference.
   */
   inline EnergyEnsemble& Simulation::energyEnsemble()
   {
      assert(energyEnsemblePtr_);
      return *energyEnsemblePtr_;
   }

   /*
   * Get the BoundaryEnsemble by reference.
   */
   inline BoundaryEnsemble& Simulation::boundaryEnsemble()
   {
      assert(boundaryEnsemblePtr_);
      return *boundaryEnsemblePtr_;
   }

   /*
   * Get the FileMaster by reference.
   */
   inline FileMaster& Simulation::fileMaster()
   {
      assert(fileMasterPtr_);
      return *fileMasterPtr_;
   }

   /*
   * Get the DiagnosticManager by reference.
   */
   inline DiagnosticManager& Simulation::diagnosticManager()
   {
      assert(diagnosticManagerPtr_);
      return *diagnosticManagerPtr_;
   }

   /*
   * Get maximum number of atom types.
   */
   inline int Simulation::nAtomType()
   {   return nAtomType_; }

   /*
   * Get maximum number of bond types.
   */
   inline int Simulation::nBondType()
   {  return nBondType_; }

   #if INTER_ANGLE
   /*
   * Get maximum number of angle types.
   */
   inline int Simulation::nAngleType()
   {  return nAngleType_; }
   #endif

   #if INTER_DIHEDRAL
   /*
   * Get maximum number of dihedral types.
   */
   inline int Simulation::nDihedralType()
   {  return nDihedralType_; }
   #endif

   #if INTER_EXTERNAL
   /*
   * Does this simulation have an external potential?
   */
   inline bool Simulation::hasExternal()
   {  return hasExternal_; }
   #endif

   /*
   * Get an AtomType descriptor by reference.
   */
   inline AtomType& Simulation::atomType(int i)
   {  return atomTypes_[i]; }

   inline MaskPolicy Simulation::maskedPairPolicy() const
   {  return maskedPairPolicy_; }

   inline bool Simulation::forceCommFlag() const
   {  return forceCommFlag_; }

}
#endif
