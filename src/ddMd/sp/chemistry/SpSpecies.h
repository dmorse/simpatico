#ifndef DDMD_SP_SPECIES_H
#define DDMD_SP_SPECIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpMolecule.h"
#include <util/containers/DArray.h>
#include <util/containers/DSArray.h>
#include <util/containers/ArrayIterator.h>

namespace DdMd
{

   class SpAtom;
   using namespace Util;

   /**
   * A set of identical molecules.
   *
   * \ingroup DdMd_Sp_Chemistry_Module
   */
   class SpSpecies {
   public:
   
     typedef ArrayIterator<SpMolecule> SpMoleculeIterator;

     /**
     * Constructor.
     */
     SpSpecies();
  
     /**
     * Set the index for this species.
     *
     * \param id species index
     */ 
     void setId(int id);
   
     /**
     * Initialize memory for this species and set sizes.
     *
     * \param nAtom number of atoms per molecule (exact)
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
     * \param atom SpAtom object to be added, must have AtomContext info.
     */ 
     void addAtom(SpAtom& atom);

     /**
     * Clear all molecules, reset to empty.
     */ 
     void clear();
  
     /**
     * Initialize an iterator over molecules.
     */ 
     void begin(SpMoleculeIterator& iterator);
  
     /**
     * Return a specific molecule by reference
     *
     * \param i molecule index
     */ 
     SpMolecule& molecule(int i);
  
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
     DArray<SpAtom*> atomPtrs_;

     // Array of molecules, each associated with a block in atomPtrs_.
     DSArray<SpMolecule> molecules_;
   
     /// SpSpecies index.
     int id_;
    
     /// Number of atoms per molecule.
     int nAtom_;
   
     /// Maximum number of molecules in this species.
     int capacity_;
   
   //friends:
   
     friend std::istream& operator >> (std::istream& in, SpSpecies& species);
     friend std::ostream& operator << (std::ostream& out, const SpSpecies& species);
   
   };

   /**
   * istream extractor (>>) for a SpSpecies.
   *
   * \param in  input stream
   * \param species  SpSpecies to be read from stream
   * \return modified  input stream
   */
   std::istream& operator>>(std::istream& in, SpSpecies &species);

   /**
   * ostream inserter (<<) for a SpSpecies.
   *
   * Format on one line with no line break:
   *
   * \param out  output stream
   * \param species  SpSpecies to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const SpSpecies &species);

   // Inline function definitions

   // Return a specific molecule. 
   inline SpMolecule& SpSpecies::molecule(int i)
   {
      return molecules_[i]; 
   }
 
   // Return integer id for this species.
   inline int SpSpecies::id() const
   {  return id_; }
 
   // Return number of atoms per molecule.
   inline int SpSpecies::nAtom() const
   {  return nAtom_; }
 
   // Number of molecules.
   inline int SpSpecies::size() const
   {  return molecules_.size(); }
 
   // Maximum number of molecules.
   inline int SpSpecies::capacity() const
   {  return capacity_; }
 
}
#endif
