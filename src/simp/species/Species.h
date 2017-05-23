#ifndef SIMP_SPECIES_H
#define SIMP_SPECIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/containers/DArray.h>      // member template
#include <util/containers/FSArray.h>     // member template
#include <simp/species/SpeciesGroup.h>   // member template parameter

namespace McMd {
   class SpeciesMutator;
}

namespace Simp
{

   using namespace Util;


   /**
   * A Species represents a set of chemically similar molecules.
   *
   * A Species has:
   *
   *    - A description of the structure of a molecule of this species.
   *
   *    - A capacity, which is the maximum allowed number of Molecules of
   *      this Species, or the number for which space should be allocated.
   *
   * The chemical structure of each molecule within a Species is defined
   * by specifying the type of each atom and constructing a list of bond,
   * angle and dihedral covalent groups. Information about covalent groups 
   * is stored internally in a set of SpeciesGroup<N> objects with N=2,3, 
   * and 4. This chemical structure is defined within the virtual protected 
   * readSpeciesParam(std::istream&) function, which is called by readParam.
   *
   * Implementations of readSpeciesParam() must read whatever 
   * parameters are required to define the chemical structure of some 
   * class of molecular structures. The default implementation of 
   * Species::readSpeciesParam() is general enough to describe any molecular 
   * structure, and simply reads all of the required information about atom 
   * types and groups from file. Subclasses of Species generally hard-code
   * some or all of the information required to describe a structure into 
   * the class definition, and require correspondingly less information to
   * to be specified in the parameter file. For example, the parameter file
   * format for a class that represents a flexible linear homopolymer might
   * require only the number of monomers per chain, a single atom type 
   * index and a single bond type index. A subclass that represents a
   * specific chemical structure (e.g., a water molecule) might not require
   * any data from the parameter file.
   *
   * A subclass of Species that is defined in the McMd namespace may also 
   * represent a "mutable" species. A mutable species is one in which each 
   * molecule can be in any of a finite number of different internal states, 
   * in which the chemical structure is slightly different in different states. 
   * A mutable species could, for example, be used implement a semi-grand 
   * ensemble for a polymer mixture, in which each polymer can be of either 
   * of two types.  A subclass that represents a mutable species must have 
   * an associated McMd::SpeciesMutator object. See documentation for the
   * Species::setMutatorPtr() function.
   *
   * \sa \ref simp_species_Species_page "parameter file format"
   * \ingroup Simp_Species_Module
   */
   class Species : public ParamComposite
   {

   public:

      #ifdef SIMP_BOND
      /// Maximum number of bonds that can be connected to one atom.
      static const int MaxBondPerAtom = 4;

      /// A SpeciesBond has the local atom ids and a type id for one bond.
      typedef SpeciesGroup<2> SpeciesBond;

      /// An array of local integer bond ids for all bonds containing one atom.
      typedef FSArray<int, MaxBondPerAtom> AtomBondIdArray;
      #endif

      #ifdef SIMP_ANGLE
      /// Maximum number of angles groups that can contain one atom.
      static const int MaxAnglePerAtom = 18;

      /// A SpeciesAngle has local atom ids and a type id for one angle.
      typedef SpeciesGroup<3> SpeciesAngle;

      /// An array of local angle ids for all angles containing one atom.
      typedef FSArray<int, MaxAnglePerAtom> AtomAngleIdArray;
      #endif

      #ifdef SIMP_DIHEDRAL
      /// Maximum number of dihedral groups that can contain one atom.
      static const int MaxDihedralPerAtom = 72;

      /// A SpeciesDihedral has local atom ids and a type id for one dihedral. 
      typedef SpeciesGroup<4> SpeciesDihedral;

      /// An array of local angle ids for all dihedrals containing one atom.
      typedef FSArray<int, MaxDihedralPerAtom> AtomDihedralIdArray;
      #endif

      /**
      * Constructor.
      */
      Species();

      /**
      * Destructor.
      */
      virtual ~Species();

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize structure for this species.
      *
      * This function reads the parameter moleculeCapacity (the maximum 
      * allowed number of molecules of this Species) and invokes the
      * virtual Species::readSpeciesParam() function to define the 
      * chemical structure of the Species.
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
      
      //@}
      /// \name Chemical Structure Accessors
      //@{
      
      /**
      * Get number of atoms per molecule for this Species.
      */
      int nAtom() const;

      /**
      * Get atom type index for a specific atom, by local atom index.
      *
      * \param iAtom local index for an atom within a molecule.
      */
      int atomTypeId(int iAtom) const;

      #ifdef SIMP_BOND
      /**
      * Get number of bonds per molecule for this Species.
      */
      int nBond() const;

      /**
      * Get a specific SpeciesBond object, by local bond index.
      *
      * \param iBond local index of a SpeciesBond within a molecule.
      */
      const SpeciesBond& speciesBond(int iBond) const;

      /**
      * Get array of ids for Bonds that contain one Atom.
      *
      * \param atomId local index for atom of interest
      */
      const AtomBondIdArray& atomBondIds(int atomId) const;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Get number of angles per molecule for this Species.
      */
      int nAngle() const;

      /**
      * Get a specific SpeciesAngle object, by local angle index.
      *
      * \param iAngle local index of a SpeciesAngle within a molecule.
      */
      const SpeciesAngle& speciesAngle(int iAngle) const;

      /**
      * Get array of ids for angles that contain one Atom.
      *
      * \param atomId local index for atom of interest
      */
      const AtomAngleIdArray& atomAngleIds(int atomId) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Get number of dihedrals per molecule for this Species.
      */
      int nDihedral() const;

      /**
      * Get a specific SpeciesDihedral object, by local angle index.
      *
      * \param iDihedral  local index of a SpeciesDihedral within a molecule.
      */
      const SpeciesDihedral& speciesDihedral(int iDihedral) const;

      /**
      * Get array of ids for dihedrals that contain one Atom.
      *
      * \param atomId  local index for atom of interest
      */
      const AtomDihedralIdArray& atomDihedralIds(int atomId) const;
      #endif

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
      McMd::SpeciesMutator& mutator();

      //@}
      /// \name Miscellaneous Accessors
      //@{

      /**
      * Get integer id of this Species.
      */
      int id() const;

      /**
      * Maximum allowed number of molecules for this Species.
      */
      int capacity() const;

      /**
      * Return true if Species is valid, or throw an Exception.
      */
      bool isValid() const;

      //@}

   protected:

      // Static constant
 
      /// Null (unknown) value for any non-negative index.
      static const int NullIndex = -1;

      // Non-static protected variables
      
      // These variables are protected, rather than private, so that they can 
      // be read and written by the (read|write)Param() functions of subclasses.

      /// Integer index for this Species.
      int id_;

      /// Number of molecules associated with the species.
      int moleculeCapacity_;

      /**
      * Number of atoms per molecule.
      */
      int nAtom_;

      /**
      * Array of atom type Ids, indexed by local atom id.
      * 
      * Element atomTypeIds_[id] is the atom type index for atom id of any
      * molecule of this species, where 0 <= id < nAtom_.
      */ 
      DArray<int> atomTypeIds_;

      #ifdef SIMP_BOND
      /**
      * Number of bonds per molecule.
      */
      int nBond_;

      /**
      * Array of SpeciesBonds for all bonds, indexed by local bond id.
      * 
      * Element speciesBonds_[id] is the SpeciesBond object for bond number id 
      * in any molecule of this Species, where 0 <= id < nBond_.
      */ 
      DArray<SpeciesBond> speciesBonds_;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Number of angles per molecule.
      */
      int nAngle_;

      /**
      * Array of SpeciesAngles for all angles, indexed by local angle id.
      * 
      * Element speciesAngles_[id] is the SpeciesAngle object for angle number id 
      * in any molecule of this Species, where 0 <= id < nAngle_.
      */ 
      DArray<SpeciesAngle> speciesAngles_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Number of dihedrals per molecule.
      */
      int nDihedral_;

      /**
      * Array of SpeciesDihedrals for all dihedrals, indexed by local dihedral id.
      * 
      * Element speciesDihedrals_[id] is the SpeciesDihedral object for dihedral
      * number id in any molecule of this Species, where 0 <= id < nAngle_.
      */ 
      DArray<SpeciesDihedral> speciesDihedrals_;
      #endif

      // Methods

      /**
      * Define chemical structure for this Species.
      *
      * This virtual function must define the structure of a molecule
      * of this species, and read any data required to do so. 
      * The default implementation Species::readSpeciesParam reads 
      * nAtom_, nBond_, nAngle_, and the elements of 
      * the arrays atomTypeIds_ and speciesBonds_, speciesAngles_.
      *
      * Re-implementations of this function by subclasses may 
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
      * This function must define the same parameter file format 
      * as readSpeciesParam.
      *
      * \param ar input parameter file stream.
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

      /**
      * Allocate chemical structure arrays.
      *
      * This function allocates the arrays that are used to define the
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

      #ifdef SIMP_BOND
      /**
      * Add a bond to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesBond object, and also adds 
      * a reference to the list of bonds that are connected to each atom.
      *
      * \param bondId   local index of bond within a molecule
      * \param atomId1  local index of 1st atom in bond
      * \param atomId2  local index of 2nd atom in bond
      * \param bondType bond type index 
      */
      void makeBond(int bondId, int atomId1, int atomId2, int bondType);
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Add an angle to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesAngle object, and also adds 
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

      #ifdef SIMP_DIHEDRAL
      /**
      * Add a dihedral to the chemical structure of a generic molecule.
      *
      * This function creates and adds a SpeciesDihedral object, and also adds 
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
      * Set a pointer to an associated McMd::SpeciesMutator for a mutable species.
      *
      * A mutable subclass of Species must have an associated SpeciesMutator 
      * object. The constructor of each such subclass should pass a pointer
      * to the SpeciesMutator to setMutatorPtr. After this function is called,
      * Species::isMutator() will return true, the associated SpeciesMutator 
      * will be accessible via the Species::mutator() function. If a mutable 
      * Species subclass is derived from both Species and SpeciesMutator 
      * (recommended), the function setMutatorPtr(this) should be invoked
      * within the subclass constructor to pass the address of the 
      * SpeciesMutator sub-object to the Species subobject.
      *
      * \param mutatorPtr pointer to an associated SpeciesMutator object
      */
      void setMutatorPtr(McMd::SpeciesMutator* mutatorPtr);

   private:

      #ifdef SIMP_BOND
      /**
      * Array of AtomBondIdArray objects for all atoms.
      * 
      * Element atomBondIdArrays_[id] is the AtomBondIdArray containing 
      * the local bond ids for bonds connected to atom with local index 
      * number i.
      */ 
      DArray<AtomBondIdArray> atomBondIdArrays_;
      #endif

      #ifdef SIMP_ANGLE
      /**
      * Array of AtomAngleIdArray objects for all atoms.
      * 
      * Element atomAngleIdArrays_[i] is the AtomAngleIdArray containing
      * the local angle ids for angles connected to atom with local index 
      * number i.
      */ 
      DArray<AtomAngleIdArray> atomAngleIdArrays_;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Array of AtomDihedralIdArray objects for all atoms.
      * 
      * Element atomDihedralIdArrays_[id] is the AtomAngleIdArray containing 
      * the local dihedral ids for dihedrals connected to atom with local 
      * index number i.
      */ 
      DArray<AtomDihedralIdArray> atomDihedralIdArrays_;
      #endif

      /**
      * Pointer to associated SpeciesMutator (if any).
      *
      * This pointer is initialized to null (0). The function isMutator()
      * returns false iff mutatorPtr_ is null. 
      */
      McMd::SpeciesMutator* mutatorPtr_;

   };

   // Inline member function definitions

   /*
   * Get integer id for this Species.
   */
   inline int Species::id() const
   { return id_; }

   /*
   * Get maximum number of molecules for this Species.
   */
   inline int Species::capacity() const
   { return moleculeCapacity_; }

   /*
   * Get number of Atoms per molecule.
   */
   inline int Species::nAtom() const
   {  return nAtom_; }

   /*
   * Get type index for atom number iAtom.
   */
   inline int Species::atomTypeId(int iAtom) const
   {  return atomTypeIds_[iAtom]; }

   #ifdef SIMP_BOND
   /*
   * Get number of Bonds per molecule
   */
   inline int Species::nBond() const
   {  return nBond_; }

   /*
   * Get a specific SpeciesBond object by local index.
   */
   inline 
   const Species::SpeciesBond& Species::speciesBond(int iBond) const
   {  return speciesBonds_[iBond]; }

   /*
   * Get array of ids for bonds that contain one atom.
   */
   inline 
   const Species::AtomBondIdArray& Species::atomBondIds(int atomId) const
   {  return atomBondIdArrays_[atomId]; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Get number of angles per molecule.
   */
   inline int Species::nAngle() const
   {  return nAngle_; }

   /*
   * Get a specific SpeciesAngle object by local index.
   */
   inline 
   const Species::SpeciesAngle& Species::speciesAngle(int iAngle) const
   {  return speciesAngles_[iAngle]; }

   /*
   * Get array of ids for angles that contain one Atom.
   */
   inline 
   const Species::AtomAngleIdArray& Species::atomAngleIds(int atomId) const
   {  return atomAngleIdArrays_[atomId]; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Get number of dihedrals per molecule
   */
   inline int Species::nDihedral() const
   {  return nDihedral_; }

   /*
   * Get a specific SpeciesDihedral object by local index.
   */
   inline 
   const Species::SpeciesDihedral& Species::speciesDihedral(int iDihedral) const
   {  return speciesDihedrals_[iDihedral]; }

   /*
   * Get array of ids for dihedrals that contain one Atom.
   */
   inline 
   const Species::AtomDihedralIdArray& Species::atomDihedralIds(int atomId) const
   {  return atomDihedralIdArrays_[atomId]; }
   #endif

   /*
   * Is this a mutable species?
   */
   inline bool Species::isMutable() const
   {  return (mutatorPtr_ != 0); }

   /*
   * Get associated mutator object.
   */
   inline McMd::SpeciesMutator& Species::mutator()
   {
      assert(mutatorPtr_);  
      return *mutatorPtr_; 
   }

} 
#endif
