/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"
#include <util/containers/RArray.h>

namespace McMd
{

   using namespace Util;

   // Define and initialize static variables
   Atom*      Atom::atoms_ = 0;
   Mask*      Atom::masks_ = 0;
   Molecule** Atom::moleculePtrs_ = 0;
   Vector*    Atom::velocities_ = 0;
   Vector*    Atom::forces_ = 0;
   bool*      Atom::isActives_ = 0;
   #ifdef MCMD_SHIFT
   IntVector* Atom::shifts_ = 0;
   #endif
   int        Atom::capacity_ = 0;

   // Static functions

   /*
   * Initialize static variables
   *
   * The only purpose of this function is to guarantee inclusion of file in static library.
   */
   void Atom::initStatic()
   {
      atoms_ = 0;
      masks_ = 0;
      moleculePtrs_ = 0;
      velocities_ = 0;
      forces_ = 0;
      isActives_ = 0;
      #ifdef MCMD_SHIFT
      shifts_ = 0;
      #endif
      capacity_ = 0;
   }

   /*
   * Allocate a static array of Atom objects.
   */
   void Atom::allocate(int capacity, RArray<Atom>& atoms)
   {
      if (capacity == 0) return;

      assert(atoms_ == 0);
      assert(masks_ == 0);
      assert(moleculePtrs_ == 0);
      assert(velocities_ == 0);
      assert(forces_ == 0);
      assert(isActives_ == 0);
      #ifdef MCMD_SHIFT
      assert(shifts_ == 0);
      #endif

      atoms_ = new Atom[capacity];
      masks_ = new Mask[capacity];
      moleculePtrs_ = new Molecule*[capacity];
      capacity_ = capacity;
      for (int i = 0; i < capacity_; ++i) {
         atoms_[i].id_ = i;
         moleculePtrs_[i] = 0;
      }
      atoms.associate(atoms_, capacity_);

      velocities_ = new Vector[capacity_];
      for (int i = 0; i < capacity_; ++i) {
         velocities_[i].zero();
      }

      forces_ = new Vector[capacity_];
      for (int i = 0; i < capacity_; ++i) {
         forces_[i].zero();
      }

      isActives_ = new bool[capacity_];
      for (int i = 0; i < capacity_; ++i) {
         isActives_[i] = true;
      }

      #ifdef MCMD_SHIFT
      shifts_ = new IntVector[capacity_];
      for (int i = 0; i < capacity_; ++i) {
         shifts_[i].zero();
      }
      #endif

   }

   /*
   * Deallocate all atoms.
   */
   void Atom::deallocate()
   {
      if (atoms_) {
         delete [] atoms_;
         atoms_ = 0;
      }
      if (masks_) {
         delete [] masks_;
         masks_ = 0;
      }
      if (moleculePtrs_) {
         delete [] moleculePtrs_;
         moleculePtrs_ = 0;
      }
      if (velocities_) {
         delete [] velocities_;
         velocities_ = 0;
      }
      if (forces_) {
         delete [] forces_;
         forces_ = 0;
      }
      if (isActives_) {
         delete [] isActives_;
         isActives_ = 0;
      }
      #ifdef MCMD_SHIFT
      if (shifts_) {
         delete [] shifts_;
         shifts_ = 0;
      }
      #endif
      capacity_ = 0;
   }

   // Nonstatic member functions

   /*
   * Constructor.
   */
   Atom::Atom() :
     typeId_(NullIndex),
     id_(NullIndex)
   {}

   /*
   * Set pointer to parent molecule.
   */
   void Atom::setMolecule(Molecule &molecule)
   {  moleculePtrs_[id_] = &molecule; }

}
