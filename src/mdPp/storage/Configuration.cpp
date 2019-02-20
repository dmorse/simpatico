/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Configuration.h"

namespace MdPp 
{

   /*
   * Constructor.
   */
   Configuration::Configuration()
    : 
      nSpecies_(0),
      atomCapacity_(0)
      #ifdef SIMP_BOND
      , bondCapacity_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleCapacity_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralCapacity_(0)
      #endif
      #ifdef SIMP_IMPROPERT
      , improperCapacity_(0)
      #endif
   {  setClassName("Configuration"); }

   /*
   * Destructor.
   */
   Configuration::~Configuration()
   {}

   /*
   * Open, read and close a parameter file.
   */
   void Configuration::readParam(const char* filename)
   {
      std::ifstream in;
      in.open(filename);
      readParam(in);
      in.close();
   }

   /*
   * Read parameters from file.
   */
   void Configuration::readParameters(std::istream& in)
   {
      #if 0
      // Optionally read species info
      nSpecies_ = 0; // default value
      readOptional<int>(in, "nSpecies", nSpecies_);
      if (nSpecies_ > 0) {
         species_.allocate(nSpecies_);
         for (int i = 0; i < nSpecies_; ++i) {
            species_[i].setId(i);
         }
      }
      #endif

      // Optionally read atom capacity
      atomCapacity_ = 0; // default value
      readOptional<int>(in, "atomCapacity", atomCapacity_); 
      if (atomCapacity_ > 0) {
         atoms_.allocate(atomCapacity_);
      }

      #ifdef SIMP_BOND
      bondCapacity_ = 0; // default value
      readOptional<int>(in, "bondCapacity", bondCapacity_); 
      // If bondCapacity is absent, it is set to zero by default
      if (bondCapacity_ > 0) {
         bonds_.allocate(bondCapacity_);
      }
      #endif

      #ifdef SIMP_ANGLE
      angleCapacity_ = 0; // default value
      readOptional<int>(in, "angleCapacity", angleCapacity_); 
      if (angleCapacity_ > 0) {
         angles_.allocate(angleCapacity_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      dihedralCapacity_ = 0; // default value
      readOptional<int>(in, "dihedralCapacity", dihedralCapacity_); 
      if (dihedralCapacity_ > 0) {
         dihedrals_.allocate(dihedralCapacity_);
      }
      #endif

      #ifdef SIMP_IMPROPER
      improperCapacity_ = 0; // default value
      readOptional<int>(in, "improperCapacity", improperCapacity_); 
      if (improperCapacity_ > 0) {
         impropers_.allocate(improperCapacity_);
      }
      #endif

   }

   void Configuration::setNSpecies(int nSpecies)
   {
      UTIL_CHECK(nSpecies > 0);
      nSpecies_ = nSpecies;
      species_.allocate(nSpecies_);
      for (int i = 0; i < nSpecies_; ++i) {
         species_[i].setId(i);
      }
   }

   /*
   * Remove all atoms and groups - set to empty state.
   */
   void Configuration::clear()
   {
      UTIL_CHECK(nSpecies_ = species_.capacity());
      if (nSpecies_ > 0) {
         for (int i = 0; i < nSpecies_; ++i) {
            species(i).clear();
         }
      }
      if (atoms_.capacity() > 0) {
         atoms_.clear();
      }
      #ifdef SIMP_BOND
      if (bonds_.capacity() > 0) {
         bonds_.clear();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angles_.capacity() > 0) {
         angles_.clear();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedrals_.capacity() > 0) {
         dihedrals_.clear();
      }
      #endif
      #ifdef SIMP_IMPROPER
      if (impropers_.capacity() > 0) {
         impropers_.clear();
      }
      #endif
   }

}
