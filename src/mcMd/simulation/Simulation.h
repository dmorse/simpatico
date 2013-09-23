#ifndef MCMD_SIMULATION_H
#define MCMD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <util/param/ParamComposite.h> // base class
#include <mcMd/chemistry/Atom.h>       // member container template argument
#include <mcMd/chemistry/Molecule.h>   // member container template argument
#include <mcMd/chemistry/Bond.h>       // typedef
#ifdef INTER_ANGLE
#include <mcMd/chemistry/Angle.h>      // typedef
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/chemistry/Dihedral.h>   // typedef
#endif
#include <mcMd/chemistry/MaskPolicy.h> // member 
#ifndef MCMD_NOATOMTYPES
#include <mcMd/chemistry/AtomType.h>   // member container template parameter
#endif

#include <util/misc/FileMaster.h>      // member
#include <util/random/Random.h>        // member
#include <util/containers/RArray.h>    // member container for Atoms
#include <util/containers/DArray.h>    // member containers (Molecules, Bonds, ...)

namespace Util 
{
   template <typename T> class ArraySet;
   template <typename T> class Factory;
   class Vector;
}

namespace McMd
{

   using namespace Util;

   class Species;
   class SpeciesManager;

   class Diagnostic;
   class DiagnosticManager;

   /**
   * The main object in a simulation, which coordinates others.
   *
   * A Simulation object is the main object in a simulation. Which is 
   * instantiated in the main program, and coordinates creation and actions 
   * of other objects. Simulation is a base class, which has subclasses 
   * McSimulation and MdSimulation designed for MC and MD simulations. 
   *
   * A Simulation has the following publicly accessible members:
   *
   *  - An array of AtomType descriptor objects.
   * 
   *  - A SpeciesManager, which has one or more Species objects.
   *
   *  - A DiagnosticManager, which has one or more Diagnostic objects.
   *
   *  - A FileMaster to manage associated input and output files.
   *
   *  - A random number generator object (an instances of Util::Random).
   *
   *  - An MPI communicator (if compiled with UTIL_MPI defined).
   *
   * A Simulation also contains several global data structures.
   *
   *  - An array of all allocated Atom objects.
   *
   *  - An array of all allocated Molecule objects, for all Species.
   *
   *  - Arrays of all Bond, Angle, and Dihedral Group objects.
   *
   * Each Molecule is associated with a block of Atom objects, and blocks 
   * of Group objects (Bond, Angle, and Dihedral). These associations are
   * created during initialization, and are permanent. 
   *
   * Each subclass of Simulation also has one or more System objects.
   * An MdSimulation has one MdSystem. An McSimulation has one McSystem.
   * A System represents a set of molecules surrounded by a Boundary.
   * A System holds a set of pointers to Molecule objects in the array
   * of Molecules owned by parent Simulation. 
   *
   * \ingroup McMd_Simulation_Module
   */
   class Simulation : public ParamComposite
   {

   public:

      #ifdef UTIL_MPI
      /**
      * Constructor.
      */
      Simulation(MPI::Intracomm& communicator);
      #endif

      /**
      * Constructor.
      */
      Simulation();

      /**
      * Destructor.
      */
      virtual ~Simulation();

      /**
      * Read parameter file and initialize.
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /// \name Initialization 
      //@{

      #ifdef UTIL_MPI
      /**
      * Set MPI job to read one parameter file and one command file.
      *
      * Call this method before readParam() to allow a single parameter
      * and command file to be used in multi-processor free energy 
      * perturbation simulations, and to read the parameter file from
      * std in. It calls Util::ParamComposite:setIoCommunicator()
      * and Util::FileMaster::setParamFileStdIn().
      *
      * \param communicator MPI communicator used for parameter file.
      */
      virtual void setIoCommunicator(MPI::Intracomm& communicator);

      /**
      * Set MPI job to read one parameter file and one command file.
      *
      * Equivalent to Simulation::setIoCommunicator(communicator()).
      */
      void setIoCommunicator();
      #endif

      /**
      * Allocate and initialize a molecule set for one Species.
      * 
      * This method is called during initialization by the readParam
      * method of an associated System.
      *
      * \param set       molecule set for one Species in a System.
      * \param speciesId integer index of the relevant Species.
      */
      void 
      allocateMoleculeSet(Util::ArraySet<Molecule> &set, int speciesId) const;

      //@}
      /// \name Read-write accessors (return by non-const reference)
      //@{

      /**
      * Get the random number generator by reference.
      */
      Random& random();

      /**
      * Get a specific Species by reference.
      * 
      * \param i integer index of desired Species
      * \return reference to species with index i
      */ 
      Species& species(int i);

      /**
      * Return the Diagnostic factory by reference.
      */
      Factory<Diagnostic>& diagnosticFactory();

      /**
      * Return the Species Factory by reference.
      */
      Factory<Species>& speciesFactory();

      /**
      * Get the FileMaster object. 
      */
      FileMaster& fileMaster();

      #ifdef UTIL_MPI
      /**
      * Get the MPI communicator by reference
      */
      MPI::Intracomm& communicator();
      #endif

      //@}
      /// \name Read-only accessors (return by value or const reference)
      //@{

      /**
      * Get the number of atom types.
      */
      int nAtomType() const;

      /**
      * Get the number of bond types.
      */
      int nBondType() const;

      #ifdef INTER_ANGLE
      /**
      * Get the number of angle types.
      */
      int nAngleType() const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get the number of dihedral types.
      */
      int nDihedralType() const;
      #endif

      #ifdef MCMD_LINK
      /**
      * Get the number of link types.
      */
      int nLinkType() const;
      #endif

      #ifdef INTER_EXTERNAL
      /**
      * Does an external potential exist?
      */
      int hasExternal() const;
      #endif

      #ifdef INTER_TETHER
      /**
      * Does a tether potential exist?
      */
      int hasTether() const;
      #endif

      /**
      * Get the number of Systems in this Simulation.
      *
      * This will return 1 for an McSimulation or MdSimulation.
      */
      int nSystem() const;

      /**
      * Get the number of Species in this Simulation.
      */
      int nSpecies() const;
     
      /**
      * Get the total number of Molecules allocated.
      */
      int moleculeCapacity() const;

      /**
      * Get the total number of Atoms allocated.
      */
      int atomCapacity() const;

      /**
      * Get the total number of Bonds allocated.
      */
      int bondCapacity() const;

      #ifdef INTER_ANGLE
      /**
      * Get the total number of Angles allocated.
      */
      int angleCapacity() const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get the total number of Dihedrals allocated.
      */
      int dihedralCapacity() const;
      #endif

      /**
      * Return the value of the mask policy (MaskNone or MaskBonded).
      */
      MaskPolicy maskedPairPolicy() const;

      /**
      * Get a specific Species by const reference.
      * 
      * \param i integer index of desired Species
      */ 
      const Species& species(int i) const;

      #ifndef MCMD_NOATOMTYPES
      /**
      * Get a single AtomType object by const reference.
      *
      * \param i integer index of desired AtomType
      */
      const AtomType& atomType(int i) const;

      /**
      * Get a const Array of all AtomType objects.
      */
      const Array<AtomType>& atomTypes() const;
      #endif

      /**
      * Get value of step index for main MC or MD loop.
      */
      int iStep() const;

      /**
      * Return true if Simulation is valid, or throw an Exception.
      */
      virtual bool isValid() const;

      //@}

   protected:

      /**
      * Number of Systems of interacting molecules (> 1 in Gibbs ensemble).
      *
      * Protected so it can be read from file and modified by a Gibbs subclass.
      * Note that nSystem_ is initialized to 1 in the Simulation constructor.
      */
      int  nSystem_;

      /**
      * Step index for main MC or MD loop.
      */
      int  iStep_;

      /**
      * Set the associated DiagnosticManager.
      *
      * This is used in the constructor for each subclass of Simulation
      * (e.g., in McSimulation and MdSimulation) by register an instance 
      * of an appropriate subclass of DiagnosticManager (e.g, either a
      * McDiagnosticManager or a MdDiagnosticManager). 
      */
      void setDiagnosticManager(DiagnosticManager* ptr);

      /**
      * Get the associated DiagnosticManager by reference.
      */
      DiagnosticManager& diagnosticManager();

   private:

      /// \name Composite members
      //@{
      
      /// Random number generator.
      Random random_;

      #ifndef MCMD_NOATOMTYPES
      /**
      * Array of AtomType objects for all types in simulation.
      *
      * Each AtomType stores the name and mass for a type of Atom.
      */
      DArray<AtomType> atomTypes_;
      #endif

      /** 
      * Array of all Molecule objects, for all Species.
      */
      DArray<Molecule> molecules_;

      /**
      * Array of all Atom objects.
      *
      * The atoms_ Array contains all Atom objects available in a simulation.
      * The Atoms associated with one Molecule form a sequential block. Blocks
      * of Atoms associated with Molecules of the same Species are also stored
      * sequentially within a larger block of Atoms allocated for that Species.
      *
      * This RArray is an alias (i.e., a shallow copy) of a DArray that is a
      * private static member of the Atom class.
      */
      RArray<Atom>     atoms_;

      /**
      * Array of all Bond objects.
      *
      * The organization of bonds_ is closely anologous to of atoms_: The 
      * Bonds associated with a Molecule are stored in a contiguous block, 
      * and blocks associated with molecules are of the same Species are 
      * stored sequentially within a larger block.
      */
      DArray<Bond>     bonds_;

      #ifdef INTER_ANGLE
      /**
      * Array of all Angle objects.
      *
      * The organization of angles_ is closely anologous to that of bonds_: The 
      * Angles associated with a Molecule are stored in a contiguous block, 
      * and blocks associated with molecules are of the same Species are 
      * stored sequentially within a larger block.
      */
      DArray<Angle>     angles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Array of all Dihedral objects.
      *
      * The organization of dihedrals_ is closely anologous to that of angles_:
      * The Dihedrals associated with a Molecule are stored in a contiguous
      * block, and blocks associated with molecules are of the same Species are
      * stored sequentially within a larger block.
      */
      DArray<Dihedral>   dihedrals_;
      #endif

      /**
      * Array containing indices to the first Molecule of each species.
      *
      * Element firstAtomIds[i] is an integer index for the first Molecule of the
      * block of the molecules_ Array associated with species number i.
      */
      DArray<int>      firstMoleculeIds_;

      /**
      * Array containing indices to the first Atom of each species.
      *
      * Element firstAtomIds[i] is an integer index for the first Atom of the
      * block of the atoms_ Array associated with species number i.
      */
      DArray<int>      firstAtomIds_;

      /**
      * Array containing indices to the first Bond of each species.
      *
      * Element firstBondIds[i] is an integer index for the first Bond of the
      * block of the bonds_ Array associated with species number i.
      */
      DArray<int>      firstBondIds_;

      #ifdef INTER_ANGLE
      /**
      * Array containing indices to the first Angle of each species.
      *
      * Element firstAngleIds[i] is an integer index for the first Angle of the
      * block of the angles_ Array associated with species number i.
      */
      DArray<int>      firstAngleIds_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Array containing indices to the first Dihedral of each species.
      *
      * Element firstDihedralIds[i] is an integer index for the first Dihedral of the
      * block of the dihedrals_ Array associated with species number i.
      */
      DArray<int>      firstDihedralIds_;
      #endif

      /**
      * Object for opening associated input and output files.
      */
      FileMaster       fileMaster_;

      #ifdef UTIL_MPI
      /**
      * Stream for log file output (serial jobs use std::cout).
      */
      std::ofstream logFile_;
      #endif

      //@}
      /// \name Pointer members
      //@{

      /**
      * Manager for molecular Species.
      *
      * An instance of SpeciesManager is created in the constructor.
      */
      SpeciesManager* speciesManagerPtr_;

      /**
      * Manager for data analysis and output classes.
      *
      * An instance of a default subclass of DiagnosticManager must be 
      * instantiated by each subclass of Simulation.
      */
      DiagnosticManager* diagnosticManagerPtr_;

      #ifdef UTIL_MPI
      /// Pointer to the simulation communicator.
      MPI::Intracomm* communicatorPtr_;
      #endif
 
      //@}
      /// \name Integer and Bool members
      //@{

      /// Number of atom types.
      int  nAtomType_;

      /// Number of bond types.
      int  nBondType_;

      #ifdef INTER_ANGLE
      /// Number of angle types.
      int  nAngleType_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Number of dihedral types.
      int  nDihedralType_;
      #endif

      #ifdef MCMD_LINK
      /// Number of link types.
      int  nLinkType_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Does an external potential exist? (0 false or 1 true)
      int  hasExternal_;
      #endif

      #ifdef INTER_TETHER
      /// Does a tether potential exist? (0 false or 1 true)
      int  hasTether_;
      #endif

      /**
      * Number of molecules allocated.
      *
      * The number of Molecule objects allocated in the molecules_ Array, 
      * for all Species in all Systems. This should be equal to the sum 
      * of values of capacity() for each Species.
      */
      int  moleculeCapacity_;

      /**
      * Number of atoms allocated.
      *
      * The number of Atom objects allocated in the atoms_ Array, for all 
      * Species in all Systems.  The atomCapacity should be equal to the
      * sum of values of capacity()*nAtom() for each Species.
      */
      int  atomCapacity_;

      /**
      * Number of bonds allocated.
      *
      * The number of Bond objects allocated in the DArray bonds_ , for all
      * Species in all Systems.
      */
      int  bondCapacity_;

      #ifdef INTER_ANGLE
      /**
      * Number of angles allocated.
      *
      * The number of Angle objects allocated in the DArray angles_, for all
      * Species in all Systems.
      */
      int  angleCapacity_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Number of dihedrals allocated.
      *
      * The number of Dihedral objects allocated in the DArray dihedrals_, for all
      * Species in all Systems.
      */
      int  dihedralCapacity_;
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

      //@}

      /**
      * Initialize all private data structures.
      */
      void initialize();

      /**
      * Initialize all Molecule and Atom objects for one Species.
      *
      * \param speciesId integer Id of the Species.
      */
      void initializeSpecies(int speciesId);
   
      /**
      * Initialize all Bond objects in all Molecules of one Species.
      *
      * \param speciesId integer Id of the Species.
      */
      void initializeSpeciesBonds(int speciesId);

      #ifdef INTER_ANGLE
      /**
      * Initialize all Angle objects in all Molecules of one Species.
      *
      * \param speciesId integer Id of the Species.
      */
      void initializeSpeciesAngles(int speciesId);
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Initialize all Dihedral objects in all Molecules of one Species.
      *
      * \param speciesId integer Id of the Species.
      */
      void initializeSpeciesDihedrals(int speciesId);
      #endif

   }; // end class Simulation


   // public inline accessor methods

   inline int Simulation::nAtomType() const
   {  return nAtomType_; }

   inline int Simulation::nBondType() const
   {  return nBondType_; }

   #ifdef INTER_ANGLE
   inline int Simulation::nAngleType() const
   {  return nAngleType_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline int Simulation::nDihedralType() const
   {  return nDihedralType_; }
   #endif

   #ifdef MCMD_LINK
   inline int Simulation::nLinkType() const
   {  return nLinkType_; }
   #endif

   #ifdef INTER_EXTERNAL
   inline int Simulation::hasExternal() const
   {  return hasExternal_; }
   #endif

   #ifdef INTER_TETHER
   inline int Simulation::hasTether() const
   {  return hasTether_; }
   #endif

   inline int Simulation::nSystem() const
   {  return nSystem_; }

   inline int Simulation::moleculeCapacity() const
   {  return moleculeCapacity_; }

   inline int Simulation::atomCapacity() const
   {  return atomCapacity_; }

   inline int Simulation::bondCapacity() const
   {  return bondCapacity_; }

   #ifdef INTER_ANGLE
   inline int Simulation::angleCapacity() const
   {  return angleCapacity_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline int Simulation::dihedralCapacity() const
   {  return dihedralCapacity_; }
   #endif

   inline MaskPolicy Simulation::maskedPairPolicy() const
   {  return maskedPairPolicy_; }

   inline Random& Simulation::random()
   {  return random_; }

   #ifndef MCMD_NOATOMTYPES
   inline const AtomType& Simulation::atomType(int i) const
   {  return atomTypes_[i]; }
   #endif

   #ifdef UTIL_MPI
   inline MPI::Intracomm& Simulation::communicator()
   {
      assert(communicatorPtr_);  
      return *communicatorPtr_; 
   }
   #endif

   // Protected inline methods
   
   inline FileMaster& Simulation::fileMaster()
   {  return fileMaster_; }

   inline void Simulation::setDiagnosticManager(DiagnosticManager* ptr)
   {  diagnosticManagerPtr_ = ptr; }

   inline DiagnosticManager& Simulation::diagnosticManager()
   {  
      assert(diagnosticManagerPtr_);
      return *diagnosticManagerPtr_; 
   }

   inline int Simulation::iStep() const
   { return iStep_; }

}
#endif
