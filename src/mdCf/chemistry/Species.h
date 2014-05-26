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
#include <util/containers/DSArray.h>
#include <util/containers/ArrayIterator.h>

namespace MdCf
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
     * \param capacity number of atoms per molecule (exact)
     * \param capacity maximum number of molecules in species
     */ 
     void initialize(int nAtom, int capacity);
  
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
     * Return number of atoms per molecule.
     */
     int nAtom() const;
  
     /** 
     * Return number of molecules in this species (=maximum id + 1)
     */
     int size() const;
   
     /** 
     * Return maximum number of molecules in this species.
     */
     int capacity() const;
   
     bool isValid() const;
   
   private:
  
     // Array of pointers to atoms, ordered by molecule 
     DArray<Atom*> atomPtrs_;

     // Array of molecules, each associated with a block in atomPtrs_.
     DSArray<Molecule> molecules_;
   
     /// Species index.
     int id_;
    
     /// Number of atoms per molecule.
     int nAtom_;
   
     /// Maximum number of molecules in this species.
     int capacity_;
   
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
      return molecules_[i]; 
   }
 
   // Return integer id for this species.
   inline int Species::id() const
   {  return id_; }
 
   // Return number of atoms per molecule.
   inline int Species::nAtom() const
   {  return nAtom_; }
 
   // Number of molecules.
   inline int Species::size() const
   {  return molecules_.size(); }
 
   // Maximum number of molecules.
   inline int Species::capacity() const
   {  return capacity_; }
 
}
#endif
