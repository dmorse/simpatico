/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Configuration.h"

namespace Tools 
{

   /*
   * Constructor.
   */
   Configuration::Configuration()
    : atomCapacity_(0)
      #ifdef SIMP_BOND
      , bondCapacity_(0)
      #endif
      #ifdef SIMP_ANGLE
      , angleCapacity_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralCapacity_(0)
      , improperCapacity_(0)
      #endif
      , nSpecies_(0)
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
      read<int>(in, "atomCapacity", atomCapacity_); 

      atoms_.allocate(atomCapacity_);

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

      improperCapacity_ = 0; // default value
      readOptional<int>(in, "improperCapacity", improperCapacity_); 
      if (improperCapacity_ > 0) {
         impropers_.allocate(improperCapacity_);
      }
      #endif

      // Optionally read species info
      nSpecies_ = 0; // default value
      readOptional<int>(in, "nSpecies", nSpecies_);
      if (nSpecies_ > 0) {
         species_.allocate(nSpecies_);
         for (int i = 0; i < nSpecies_; ++i) {
            species_[i].setId(i);
         }
         readDArray<Species>(in, "species", species_, nSpecies_);
      }

   }

   /*
   * Remove all atoms and groups - set to empty state.
   */
   void Configuration::clear()
   {
      atoms_.clear();
      #ifdef SIMP_BOND
      if (bondCapacity_ > 0) {
         bonds_.clear();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleCapacity_ > 0) {
         angles_.clear();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralCapacity_ > 0) {
         dihedrals_.clear();
      }
      if (improperCapacity_ > 0) {
         impropers_.clear();
      }
      #endif
   }

}
