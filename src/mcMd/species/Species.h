#ifndef MCMD_SPECIES_H
#define MCMD_SPECIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <util/containers/ArrayStack.h>    // member template
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

#include <mcMd/chemistry/Molecule.h>       // member template parameter
#include <mcMd/chemistry/SpeciesGroup.h>   // member template parameter
#include <mcMd/chemistry/Bond.h>           // typedef
#ifdef INTER_ANGLE
#include <mcMd/chemistry/Angle.h>          // typedef
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/chemistry/Dihedral.h>       // typedef
#endif
#include <util/boundary/Boundary.h>        // typedef

namespace McMd
{

   using namespace Util;

   class Atom;
   class SpeciesMutator;
   class BondPotential;

   /**
   * A Species represents a set of chemically similar molecules.
   *
   * A Species has:
   *
   *    - A description of the structure of a molecule of this species.
   *
   *    - A capacity, which is the number of Molecules of this Species for
   *      which memory must be allocated.
   *
   *    - A reservoir, which holds pointers to unused Molecule objects.
   *
   * The chemical structure is initialized in the readParam method and may be
   * queried by several accessor methods. The capacity specifies the number 
   * of Molecules of this of this species for which space should be allocated 
   * by the parent Simulation.  The reservoir is a stack that is used to hold
   * pointers to Molecules of this Species that are not currently owned by 
   * any System.
   *
   * The chemical structure of a generic molecule within a Species is defined
   * by specifying the type of each atom and constructing a list of bond,
   * angle and dihedral covalent groups. Information about covalent groups 
   * is stored internally in a set of SpeciesGroup<N> objects with N=2,3, 
   * and 4. This chemical structure is defined within the virtual protected 
   * readSpeciesParam(std::istream&) method, which is called by readParam.
   *
   * Implementations of readSpeciesParam() may read any parameters that are 
   * required to define the chemical structure of an particular Species.
   * The default implementation of Species::readSpeciesParam() simply reads 
   * all of the required information from file. This is thus sufficiently 
   * flexible to describe any molecular species. Subclasses of Species may 
   * instead hard-code some or all of the information required to describe
   * a molecular structure into the class definition, and thus reduce the
   * amount of information that must be read from file. For example, if a 
   * subclass of Species represents a unique chemical structure, such as
   * water or some other specific molecule, then the entire structure can 
   * be hard-coded into the method implementtion. In other cases, in which
   * a subclass is used to describe a category of molecular species, a more
   * limited amount of information must be read from file in order to
   * specify a unique structure. For example, a class that represents a 
   * linear homopolymer might read information required describe the 
   * monomer and to specify the number of monomers per chain.
   *
   * A subclass of Species may also be usd to represent a "psuedo species".
   * A pseudo species is a species in which each molecule can be in any
   * of a finite number of different internal states, in which the chemical
   * structure is slightly different in different states. A pseudo-molecule
   * could, for example, be used implement a semi-grand ensemble for a
   * polymer mixture, in which each polymer can be of either of two types.
   * The readMoleculeState() and writeMoleculeState() methods are provided
   * to allow information about the internal state of each molecules in
   * such an ensemble to be read from and written to a configuration file.
   * The default implementations of these methods do nothing.
   *
   * \ingroup McMd_Species_Module
   */
   class Species : public ParamComposite
   {

   public:

      /// Maximum number of bonds that can be connected to one atom.
      static const int MaxBondPerAtom    = 4;

      /// Maximum number of angles groups that can contain one atom.
      static const int MaxAnglePerAtom   = 18;

      /// Maximum number of dihedral groups that can contain one atom.
      static const int MaxDihedralPerAtom = 72;

      /// A SpeciesBond has the local atom ids and a type id for one bond.
      typedef SpeciesGroup<2>               SpeciesBond;

      /// Array of pointers to Bond objects that contain a specific Atom.
      typedef FSArray<const Bond*, MaxBondPerAtom> AtomBondArray;

      #ifdef INTER_ANGLE
      /// A SpeciesAngle has the local atom ids and a type id for an angle.
      typedef SpeciesGroup<3>               SpeciesAngle;

      /// Array of pointers to Angle objects that contain a specific Atom.
      typedef FSArray<const Angle*, MaxAnglePerAtom> AtomAngleArray;
      #endif

      #ifdef INTER_DIHEDRAL
      /// A SpeciesDihedral has the local atom ids and a type id for an angle.
      typedef SpeciesGroup<4>               SpeciesDihedral;

      /// Array of pointers to Angle objects that contain a specific Atom.
      typedef FSArray<const Dihedral*, MaxDihedralPerAtom> AtomDihedralArray;
      #endif

      /// A stack of unused Molecule objects available for this Species.
      typedef ArrayStack<Molecule>          Reservoir;

      /**
      * Constructor.
      */
      Species();

      /**
      * Destructor.
      */
      virtual ~Species();

      /**
      * Read parameters and allocate molecules for this species.
      *
      * This method reads the parameter moleculeCapacity (the number 
      * of Molecule objects allocated for this Species) and invokes 
      * virtual Species::readSpeciesParam() method to define the chemical
      * structure of the Species, and read any parameters required to do 
      * so. The method also allocates memory for an associated reservoir.
      *
      * \param in input stream.
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

      /**
      * Set integer id for this Species.
      *
      * \param id integer id
      */
      void setId(int id);
      
      /// \name Chemical Structure Accessors
      //@{
      
      /**
      * Get number of atoms per molecule for this Species.
      */
      int nAtom() const;

      /**
      * Get number of bonds per molecule for this Species.
      */
      int nBond() const;

      #ifdef INTER_ANGLE
      /**
      * Get number of angles per molecule for this Species.
      */
      int nAngle() const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get number of dihedrals per molecule for this Species.
      */
      int nDihedral() const;
      #endif

      /**
      * Get atom type index for a specific atom, by local atom index.
      *
      * \param iAtom local index for an atom within a molecule.
      */
      int atomTypeId(int iAtom) const;

      /**
      * Get a specific SpeciesBond object, by local bond index.
      *
      * \param iBond local index of a SpeciesBond within a molecule.
      */
      const SpeciesBond& speciesBond(int iBond) const;

      /**
      * Get array of pointers to Bonds that contain an Atom.
      *
      * \param atom  atom of interest
      * \param bonds array of pointers to Bonds containing the atom
      */
      void getAtomBonds(const Atom& atom, AtomBondArray& bonds) const;

      #ifdef INTER_ANGLE
      /**
      * Get a specific SpeciesAngle object, by local angle index.
      *
      * \param iAngle local index of a SpeciesAngle within a molecule.
      */
      const SpeciesAngle& speciesAngle(int iAngle) const;

      /**
      * Get array of pointers to Angles that contain an Atom.
      *
      * \param atom   atom of interest.
      * \param angles array of pointers to Angles containing the atom.
      */
      void getAtomAngles(const Atom& atom, AtomAngleArray& angles) const;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get a specific SpeciesDihedral object, by local angle index.
      *
      * \param iDihedral  local index of a SpeciesDihedral within a molecule.
      */
      const SpeciesDihedral& speciesDihedral(int iDihedral) const;

      /**
      * Get array of pointers to Dihedrals that contain an Atom.
      *
      * \param atom     atom of interest.
      * \param dihedrals array of pointers to Dihedrals containing the atom.
      */
      void getAtomDihedrals(const Atom& atom, AtomDihedralArray& dihedrals) const;
      #endif

      /**
      * Generate random molecules
      *
      * \param nMolecule number of molecules to genearte
      * \param exclusionRadius array of exclusion radii for every atom type
      * \param system the System
      * \param bondPotentialPtr the bond potential
      * \param boundary the boundary to generate atoms in
      */
      virtual void generateMolecules(int nMolecule,
         DArray<double> exclusionRadius, System &system,
         BondPotential *bondPotentialPtr, const Boundary &boundary);

      //@}
      /// \name Miscellaneous Accessors
      //@{

      /**
      * Get integer id of this Species.
      */
      int id() const;

      /**
      * Get total number of Molecule objects allocated for this Species.
      */
      int capacity() const;

      /**
      * Get the reservoir for this Species by reference.
      * 
      * The reservoir is an ArrayStack<Molecule> containers of pointers to
      * all unused Molecule objects of this species, i.e., all Molecules
      * that are not owned by any System.
      */
      Reservoir& reservoir();

      /**
      * Return true if Species is valid, or throw an Exception.
      */
      bool isValid() const;

      //@}
      /// \name Interface for mutable species.
      //@{
      
      /**
      * Is this a mutable Species?
      */
      bool isMutable() const;

      /**
      * Return the species mutator object by reference.
      */
      SpeciesMutator& mutator();

      //@}
      
   protected:

      // Typedefs

      /// An array of local integer bond ids for all bonds containing one atom.
      typedef FSArray<int, MaxBondPerAtom>    AtomBondIdArray;

      #ifdef INTER_ANGLE
      /// An array of local integer angle ids for all angles containing one atom.
      typedef FSArray<int, MaxAnglePerAtom>   AtomAngleIdArray;
      #endif

      #ifdef INTER_DIHEDRAL
      /// An array of local integer angle ids for all dihedrals containing one atom.
      typedef FSArray<int, MaxDihedralPerAtom> AtomDihedralIdArray;
      #endif

      // Static constant
     
      /// Null (unknown) value for any non-negative index.
      static const int NullIndex = -1;

      // Non-static protected variables
      
      // These variables are protected, rather than private, so that they can 
      // be read and written by the (read|write)Param() methods of subclasses.

      /**
      * Array of atom type Ids, indexed by local atom id.
      * 
      * Element atomTypeIds_[id] is the atom type index for atom id of any
      * molecule of this species, where 0 <= id < nAtom_.
      */ 
      DArray<int>           atomTypeIds_;

      /**
      * Array of SpeciesBonds for all bonds, indexed by local bond id.
      * 
      * Element speciesBonds_[id] is the SpeciesBond object for bond number id 
      * in any Molecule of this Species, where 0 <= id < nBond_.
      */ 
      DArray<SpeciesBond>   speciesBonds_;

      #ifdef INTER_ANGLE
      /**
      * Array of SpeciesAngles for all angles, indexed by local angle id.
      * 
      * Element speciesAngles_[id] is the SpeciesAngle object for angle number id 
      * in any Molecule of this Species, where 0 <= id < nAngle_.
      */ 
      DArray<SpeciesAngle>   speciesAngles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Array of SpeciesDihedrals for all dihedrals, indexed by local dihedral id.
      * 
      * Element speciesDihedrals_[id] is the SpeciesDihedral object for dihedral
      * number id in any Molecule of this Species, where 0 <= id < nAngle_.
      */ 
      DArray<SpeciesDihedral> speciesDihedrals_;
      #endif

      /**
      * Number of atoms per molecule.
      */
      int    nAtom_;

      /**
      * Number of bonds per molecule.
      */
      int    nBond_;

      #ifdef INTER_ANGLE
      /**
      * Number of angles per molecule.
      */
      int    nAngle_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Number of dihedrals per molecule.
      */
      int    nDihedral_;
      #endif

      /// Number of molecules associated with the species.
      int    moleculeCapacity_;

      /// Integer index for this Species.
      int    id_;

      // Methods

      /**
      * Define chemical structure for this Species.
      *
      * This virtual method must define the structure of a molecule
      * of this species, and read any data required to do so. 
      * The default implementation Species::readSpeciesParam reads 
      * nAtom_, nBond_, nAngle_, and the elements of 
      * the arrays atomTypeIds_ and speciesBonds_, speciesAngles_.
      *
      * Re-implementations of this method by subclasses may 
      * hard-code some or all of the information contained in these 
      * variables, and may define a more compact file format for 
      * any parameters that are read from file.
      * 
      * \param in input parameter file stream.
      */
      virtual void readSpeciesParam(std::istream &in);

      /**
      * Define chemical structure for this Species.
      *
      * Analogous to readSpeciesParam, but reads data from a
      * Serializable::IArchive rather than an input stream.
      * This method must define the same parameter file format 
      * as readSpeciesParam.
      *
      * \param ar input parameter file stream.
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

      /**
      * Allocate chemical structure arrays.
      *
      * This method allocates the arrays that are used to define the
      * chemical structure of a generic molecule, such as atomTypeIds_, 
      * speciesBonds_, atomBondIdArrays_, speciesAngles_, etc.
      * 
      * Precondition: nAtom_, nBond_, nAngles_, etc. must have nonzero
      * values on entry.
      */
      void allocate();

      /**
      * Set the type for one atom in a generic molecule of this Species.
      *
      * \param atomId   local atom id, in range 0,..., nAtom - 1;
      * \param atomType atom type index
      */
      void setAtomType(int atomId, int atomType);

      /**
      * Add a bond to the chemical structure of a generic molecule.
      *
      * This method creates and adds a SpeciesBond object, and also adds 
      * a reference to the list of bonds that are connected to each atom.
      *
      * \param bondId   local index of bond within a molecule
      * \param atomId1  local index of 1st atom in bond
      * \param atomId2  local index of 2nd atom in bond
      * \param bondType bond type index 
      */
      void makeBond(int bondId, int atomId1, int atomId2, int bondType);

      #ifdef INTER_ANGLE
      /**
      * Add an angle to the chemical structure of a generic molecule.
      *
      * This method creates and adds a SpeciesAngle object, and also adds 
      * a reference to the list of angles that are connected to each atom.
      *
      * \param angleId   local index of angle within a molecule
      * \param atomId1   local index of 1st atom in angle
      * \param atomId2   local index of 2nd atom in angle
      * \param atomId3   local index of 3rd atom in angle
      * \param angleType angle type index 
      */
      void makeAngle(int angleId, int atomId1, int atomId2, int atomId3,
                     int angleType);
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Add a dihedral to the chemical structure of a generic molecule.
      *
      * This method creates and adds a SpeciesDihedral object, and also adds 
      * a reference to the list of dihedrals that are connected to each atom.
      *
      * \param dihedralId  local index of dihedral within a molecule
      * \param atomId1     local index of 1st atom in dihedral
      * \param atomId2     local index of 2nd atom in dihedral
      * \param atomId3     local index of 3rd atom in dihedral
      * \param atomId4     local index of 4th atom in dihedral
      * \param dihedralType dihedral type index 
      */
      void makeDihedral(int dihedralId, int atomId1, int atomId2, int atomId3,
                     int atomId4, int dihedralType);
      #endif

      /**
      * Set a pointer to an associated SpeciesMutator for a mutable species.
      *
      * A mutable species must be implemented as a subclass of both Species 
      * and SpeciesMutator. The constructor of any such class must call 
      * setMutatorPtr(this) to give the address of the SpeciesMutator subobject
      * to the Species subobject. If an instance of such a class is accessed
      * through a Species* pointer, the SpeciesMutator is then be accessible 
      * via the Species::mutator() method.
      *
      * \param mutatorPtr pointer to an associated SpeciesMutator object
      */
      void setMutatorPtr(SpeciesMutator* mutatorPtr);

   private:

      /**
      * Reservoir of molecules that are not assigned to any System.
      *
      * The reservoir_ stack holds pointers to Molecule objects within the
      * molecules_ array that are not assigned to any System.
      */
      ArrayStack<Molecule>  reservoir_;

      /**
      * Array of AtomBondIdArray objects for all atoms.
      * 
      * Element atomBondArray_[id] is the AtomBondIdArray containing the
      * local bond ids for bonds connected to atom with local index id 
      * in a generic Molecule of this Species.
      */ 
      DArray<AtomBondIdArray> atomBondIdArrays_;

      #ifdef INTER_ANGLE
      /**
      * Array of AtomAngleIdArray objects for all atoms.
      * 
      * Element atomAngleArray_[id] is the AtomAngleIdArray containing the
      * local angle ids for angles connected to atom with local index id 
      * in a generic Molecule of this Species.
      */ 
      DArray<AtomAngleIdArray> atomAngleIdArrays_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Array of AtomDihedralIdArray objects for all atoms.
      * 
      * Element atomDihedralArray_[id] is the AtomAngleIdArray containing the
      * local dihedral ids for dihedrals connected to atom with local index id 
      * in a generic Molecule of this Species.
      */ 
      DArray<AtomDihedralIdArray> atomDihedralIdArrays_;
      #endif

      /// Pointer to associated Species Mutator (if any).
      SpeciesMutator* mutatorPtr_;

   //friends:
   
      friend class PseudoSpecies; 

   };

   // Inline method definitions

   /*
   * Get integer id for this Species.
   */
   inline int Species::id() const
   { return id_; }

   /*
   * Get number of Molecule objects allocated for this Species.
   */
   inline int Species::capacity() const
   { return moleculeCapacity_; }

   /**
   * Get the reservoir for this Species by reference.
   */
   inline Species::Reservoir& Species::reservoir()
   { return reservoir_; }

   /*
   * Get number of Atoms per Molecule.
   */
   inline int Species::nAtom() const
   { return nAtom_; }

   /*
   * Get number of Bonds per Molecule
   */
   inline int Species::nBond() const
   { return nBond_; }

   #ifdef INTER_ANGLE
   /*
   * Get number of angles per Molecule
   */
   inline int Species::nAngle() const
   { return nAngle_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Get number of dihedrals per Molecule
   */
   inline int Species::nDihedral() const
   { return nDihedral_; }
   #endif

   /*
   * Get type index for atom number iAtom.
   */
   inline int Species::atomTypeId(int iAtom) const
   { return atomTypeIds_[iAtom]; }

   /*
   * Get a specific SpeciesBond object by local index.
   */
   inline const Species::SpeciesBond& Species::speciesBond(int iBond) const
   { return speciesBonds_[iBond]; }

   #ifdef INTER_ANGLE
   /*
   * Get a specific SpeciesAngle object by local index.
   */
   inline const Species::SpeciesAngle& Species::speciesAngle(int iAngle) const
   { return speciesAngles_[iAngle]; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Get a specific SpeciesDihedral object by local index.
   */
   inline const Species::SpeciesDihedral& Species::speciesDihedral(int iDihedral) const
   { return speciesDihedrals_[iDihedral]; }
   #endif

   /*
   * Is this a mutable species?
   */
   inline bool Species::isMutable() const
   {  return (mutatorPtr_ != 0); }

   /*
   * Get associated mutator object.
   */
   inline SpeciesMutator& Species::mutator()
   {
      assert(mutatorPtr_);  
      return *mutatorPtr_; 
   }

} 
#endif
