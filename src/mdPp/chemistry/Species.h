#ifndef MDPP_SPECIES_H
#define MDPP_SPECIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"
#include <util/containers/DArray.h>
#include <util/containers/ArrayIterator.h>

namespace MdPp
{

   class Atom;
   using namespace Util;

   class Species {
   public:
   
     typedef ArrayIterator<Molecule> MoleculeIterator;

     /**
     * Constructor.
     */
     Species();
  
     /**
     * Set the index for this species.
     *
     * \param id species index
     */ 
     void setId(int id);
   
     /**
     * Initialize memory for this species and set sizes.
     *
     * \param capacity maximum number of molecules in species
     * \param capacity number of atoms per molecule (exact)
     */ 
     void initialize(int capacity, int nAtom);
  
     /**
     * Initialize memory for this species.
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
     void begin(MoleculeIterator& iterator);
  
     /**
     * Return a specific molecule by reference
     *
     * \param i molecule index
     */ 
     Molecule& molecule(int i);
  
     /** 
     * Return integer id for this species.
     */
     int id() const;
  
     /** 
     * Return number of molecules in this species (=maximum id + 1)
     */
     int size() const;
   
     bool isValid() const;
   
   private:
   
     DArray<Atom*> atomPtrs_;
     DArray<Molecule> molecules_;
   
     /// Species index.
     int id_;
    
     /// Number of atoms per molecule.
     int nAtom_;
   
     /// Maximum number of molecules in this species.
     int capacity_;
   
     /// Actual number of molecules = maximum molecule id + 1.
     int size_;
   
   //friends:
   
     friend std::istream& operator >> (std::istream& in, Species& species);
     friend std::ostream& operator << (std::ostream& out, const Species& species);
   
   };

   /**
   * istream extractor (>>) for a Species.
   *
   * \param in  input stream
   * \param species  Species to be read from stream
   * \return modified  input stream
   */
   std::istream& operator>>(std::istream& in, Species &species);

   /**
   * ostream inserter (<<) for a Species.
   *
   * Format on one line with no line break:
   *
   * \param out  output stream
   * \param species  Species to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Species &species);

   // Inline function definitions

   // Return a specific molecule. 
   inline Molecule& Species::molecule(int i)
   {
      assert(i >=0);
      assert(i < size_);
      return molecules_[i]; 
   }
 
   // Return integer id for this species.
   inline int Species::id() const
   {  return id_; }
 
   // Number of molecules.
   inline int Species::size() const
   {  return size_; }
 
}
#endif
