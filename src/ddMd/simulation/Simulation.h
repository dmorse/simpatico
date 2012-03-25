#ifndef DDMD_DD_SIMULATION_H
#define DDMD_DD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class
#include <ddMd/communicate/Domain.h>             // member 
#include <ddMd/communicate/Buffer.h>             // member 
#include <ddMd/communicate/Exchanger.h>          // member 
#include <ddMd/storage/AtomStorage.h>            // member 
#include <ddMd/storage/BondStorage.h>            // member 
#include <ddMd/chemistry/AtomType.h>             // member
#include <ddMd/chemistry/MaskPolicy.h>           // member
#include <ddMd/diagnostics/DiagnosticManager.h>  // member
#include <util/boundary/Boundary.h>              // member 
#include <util/random/Random.h>                  // member 
#include <util/containers/DArray.h>              // member 

namespace Util
{
   template <typename T> class Factory;
}


namespace DdMd
{

   class EnergyEnsemble;
   class BoundaryEnsemble;
   class PairPotential;
   class BondPotential;
   class Integrator;
   class ConfigIo;
   class FileMaster;
   class DiagnosticManager;

   using namespace Util;

   /**
   * Main object for a domain-decomposition MD simulation.
   *
   * A DdMd::Simulation contains and coordinates all the components of a parallel
   * MD simulation. 
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
      * Compute forces for all local atoms.
      */
      void computeForces();

      /**
      * Is exchange of atoms among processors needed?
      *
      * This method returns true if any atom in the AtomStorage has moved
      * a distance skin/2 or greater since the last snapshot (i.e., the
      * last time atoms were exchanged and pair list was rebuilt).
      *
      * Reduce-to-all operation: Must be called on all nodes and returns 
      * same result on all node.
      * 
      * \return true if exchange / reneighboring is needed, false otherwise.
      */
      bool needExchange();

      /**
      * Compute total kinetic energy.
      * 
      * Reduce operation: Must be called on all nodes.
      */
      void computeKineticEnergy();

      /**
      * Compute precomputed total kinetic energy.
      * 
      * Call only on master node. 
      *
      * \return total kinetic energy for all nodes on master.
      */
      double kineticEnergy();

      /**
      * Calculate total potential energy on all processors.
      * 
      * Reduce operation: Must be called on all nodes .
      *
      * \return total pair potential for all nodes on master, 0.0 otherwise.
      */
      void computePotentialEnergies();

      /**
      * Return precomputed total potential energy.
      *
      * Call only on master processor.
      * 
      * \return total potential energy (call only on master node).
      */
      double potentialEnergy();

      #if 0
      /**
      * Calculate total bond potential energy.
      * 
      * Reduce operation: Must be called on all nodes but returns correct
      * total value only on grid communicator master.
      *
      * \return total bond potential for all nodes on master, 0.0 otherwise.
      */
      double bondPotentialEnergy();
      #endif

      /// \name Config File IO
      //@{

      /**
      * Read configuration file on master and distribute atoms.
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

      #if 0
      /**
      * Get the configuration file reader/writer factory by reference.
      */
      Factory<ConfigIo>& configIoFactory();
      #endif

      //@}

      /// \name Potential Energy Factories and Styles
      //@{

      #ifndef DDMD_NOPAIR
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

      #ifdef DDMD_ANGLE
      /**
      * Get the associated AngleFactory by reference.
      */
      Factory<AnglePotential>& angleFactory();

      /**
      * Return angle potential style string.
      */
      std::string angleStyle() const;
      #endif

      #ifdef DDMD_DIHEDRAL
      /**
      * Get the associated Dihedral Factory by reference.
      */
      Factory<DihedralPotential>& dihedralFactory();

      /**
      * Return dihedral potential style string.
      */
      std::string dihedralStyle() const;
      #endif

      //@}
      /// \name Accessors (Miscellaneous)
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
      * Get maximum number of atom types.
      */
      int nAtomType();

      /**
      * Get maximum number of bond types.
      */
      int nBondType();

      /**
      * Get an AtomType descriptor by reference.
      */
      AtomType& atomType(int i);

      /**
      * Return the value of the mask policy (MaskNone or MaskBonded).
      */
      MaskPolicy maskedPairPolicy() const;

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

      /**
      * Return associated DiagnosticManager by reference.
      */
      DiagnosticManager& diagnosticManager();

   private:

      /// Container for all atoms and ghosts.
      AtomStorage   atomStorage_;

      /// Container for all atoms and ghosts.
      BondStorage   bondStorage_;

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

      #ifdef UTIL_MPI
      MPI::Intracomm* communicatorPtr_;
      #endif

      /// Pointer to force/energy evaluator.
      PairPotential* pairPotentialPtr_;

      /// Pointer to force/energy evaluator.
      BondPotential* bondPotentialPtr_;

      /// Pointer to MD integrator.
      Integrator*   integratorPtr_;

      /// Pointer to configuration file reader/writer.
      ConfigIo*     configIoPtr_;

      /// Pointer to an EnergyEnsemble.
      EnergyEnsemble*   energyEnsemblePtr_;

      /// Pointer to an BoundaryEnsemble.
      BoundaryEnsemble* boundaryEnsemblePtr_;

      /// Pointer to a FileMaster.
      FileMaster*         fileMasterPtr_;

      /// DiagnosticManager
      DiagnosticManager*  diagnosticManagerPtr_;

      #ifndef DDMD_NOPAIR
      /// Pointer to a PairPotential factory.
      Factory<PairPotential>*  pairFactoryPtr_;
      #endif

      /// Pointer to a Factory<BondPotential>.
      Factory<BondPotential>*  bondFactoryPtr_;

      #ifdef DDMD_ANGLE
      /// Pointer to the AnglePotential Factory.
      Factory<AnglePotential>*  angleFactoryPtr_;
      #endif

      #ifdef DDMD_DIHEDRAL
      /// Pointer to DihedralPotential Factory
      Factory<DihedralPotential>*  dihedralFactoryPtr_;
      #endif

      /// Pointer to MD integrator factory.
      Factory<Integrator>* integratorFactoryPtr_;

      #if 0
      /// Pointer to a configuration reader/writer factory.
      Factory<ConfigIo>* configIoFactoryPtr_;
      #endif

      #ifndef DDMD_NOPAIR
      /// Name of pair potential style.
      std::string pairStyle_;
      #endif

      /// Name of bond potential style.
      std::string bondStyle_;

      #ifdef DDMD_ANGLE
      /// Name of angle potential style.
      std::string angleStyle_;
      #endif

      #ifdef DDMD_DIHEDRAL
      /// Name of dihedral potential style.
      std::string dihedralStyle_;
      #endif

      /// Number of distinct atom types.
      int nAtomType_;

      /// Number of distinct bond types.
      int nBondType_;

      /**
      * Policy for suppressing pair interactions for some atom pairs.
      *
      * Allowed values of enum MaskPolicy:
      *
      *  - MaskNone:   no masked pairs
      *  - MaskBonded:  mask pair interaction between bonded atoms
      */
      MaskPolicy  maskedPairPolicy_;

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

   inline Exchanger& Simulation::exchanger()
   { return exchanger_; }

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
   * Get the FileMaster by reference.
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

   /*
   * Get an AtomType descriptor by reference.
   */
   inline AtomType& Simulation::atomType(int i)
   {  return atomTypes_[i]; }

   inline MaskPolicy Simulation::maskedPairPolicy() const
   {  return maskedPairPolicy_; }


}
#endif
