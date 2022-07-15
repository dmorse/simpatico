#ifndef MDPP_SPECIES_STORAGE_H
#define MDPP_SPECIES_STORAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <simp/species/Species.h>           //  Base class
#include <mdPp/chemistry/Molecule.h>        //  Template argument
#include <util/containers/DSArray.h>        //  Member template
#include <util/containers/DArray.h>         //  Member template
#include <util/containers/ArrayIterator.h>  //  template in interface

namespace MdPp
{

   struct Atom;
   using namespace Util;
   using namespace Simp;

   /**
   * A set of identical molecules.
   *
   * \ingroup MdPp_Storage_Module
   */
   class SpeciesStorage : public Species 
   {

   public:
  
      /**
      * Iterator over molecules in a species.
      */ 
      typedef ArrayIterator<Molecule> Iterator;

      /**
      * Constructor.
      */
      SpeciesStorage();
  
      /**
      * Destructor.
      */
      ~SpeciesStorage();
 
      /**
      * Initialize memory for this species storage and set sizes.
      */ 
      void initialize();
  
      /**
      * Add an atom to the species.
      *
      * \param atom Atom object to be added, must have AtomContext info.
      */ 
      void addAtom(Atom& atom);

      /**
      * Clear all molecules, reset to empty.
      */ 
      void clear();
  
      /**
      * Initialize an iterator over molecules.
      */ 
      void begin(Iterator& iterator);
  
      /**
      * Return a specific molecule by reference
      *
      * \param i molecule index
      */ 
      Molecule& molecule(int i);
  
      /** 
      * Return number of molecules in this species storage (=maximum id + 1)
      */
      int size() const;
  
      /**
      * Return true if valid, or throw Exception.
      */ 
      bool isValid() const;

      #if 0
      /**
      * Serialize to/from an archive.
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
      #endif
   
   private:

      /**  
      * Array of pointers to atoms, ordered by molecule.
      *
      * Capacity = (molecule capacity)*(# atoms per molecule)
      */
      DArray<Atom*> atomPtrs_;

      /**
      * Array of molecules, each associated with a block in atomPtrs_.
      *
      * Molecules are ordered by molecule id, as stored in the atoms.
      * The size of the array is always one greater than the maximum
      * id of a molecule for which one or more atoms have been added.
      */
      DSArray<Molecule> molecules_;
 
   //friends:
  
      #if 0 
      friend 
      std::istream& operator >> (std::istream& in, SpeciesStorage& species);

      friend 
      std::ostream& operator << (std::ostream& out, const SpeciesStorage& species);
      #endif
   
   };

   /**
   * istream extractor (>>) for a SpeciesStorage.
   *
   * \param in  input stream
   * \param species  SpeciesStorage to be read from stream
   * \return modified  input stream
   */
   std::istream& operator>>(std::istream& in, SpeciesStorage &species);

   /**
   * ostream inserter (<<) for a SpeciesStorage.
   *
   * Format on one line with no line break:
   *
   * \param out  output stream
   * \param species  SpeciesStorage to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const SpeciesStorage &species);

   // Inline function definitions

   /*
   * Return a specific molecule. 
   */
   inline Molecule& SpeciesStorage::molecule(int i)
   {  return molecules_[i]; }

   /* 
   * Return the number of molecules (maximum molecule id + 1).
   */
   inline int SpeciesStorage::size() const
   {  return molecules_.size(); }

   #if 0
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void SpeciesStorage::serialize(Archive& ar, const unsigned int version)
   {
      ar & nAtom_;
      ar & capacity_;
      if (Archive::is_loading()) {
         initialize();
      }
   }
   #endif

}
#endif
