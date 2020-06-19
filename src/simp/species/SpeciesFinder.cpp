/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2020, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesFinder.h"

namespace Simp
{

   using namespace Util;

   SpeciesFinder::SpeciesFinder()
    : nSpecies_(-1)
   {}

   void SpeciesFinder::allocate(int nSpecies)
   {
      UTIL_CHECK(nSpecies > 0);
      nSpecies_ = nSpecies;

      nMolecule_.allocate(nSpecies);   
      nPart_.allocate(nSpecies);   
      for (int i = 0; i < nSpecies_; ++i) {
         nMolecule_[i] = 0;
         nPart_[i] = 0;
      }

      firstMoleculeId_.allocate(nSpecies+1);   
      firstPartId_.allocate(nSpecies+1);   
      for (int i = 0; i <= nSpecies_; ++i) {
         firstMoleculeId_[i] = 0;
         firstPartId_[i] = 0;
      }

   }

   void SpeciesFinder::setSpecies(int iSpecies, int nMolecule, int nPart)
   {
      UTIL_CHECK(nSpecies_ > 0);
      UTIL_CHECK(iSpecies >= 0);
      UTIL_CHECK(iSpecies < nSpecies_);
      nMolecule_[iSpecies] = nMolecule;
      nPart_[iSpecies] = nPart;
   }

   void SpeciesFinder::initialize()
   {
      int nMolecule, nPart;
      firstMoleculeId_[0] = 0;
      firstPartId_[0] = 0;
      for (int i = 0; i < nSpecies_; ++i) {
         nMolecule = nMolecule_[i];
         nPart = nPart_[i];
         firstMoleculeId_[i+1] = firstMoleculeId_[i] + nMolecule;
         firstPartId_[i+1] = firstPartId_[i] + nMolecule*nPart;
      }
   }

   void SpeciesFinder::findMolecule(int moleculeId, 
                                    SpeciesFinder::Molecule& context)
   {
      UTIL_CHECK(nSpecies_ > 0);
      UTIL_CHECK(moleculeId >= 0);
      UTIL_CHECK(moleculeId < firstMoleculeId_[nSpecies_]);
      bool found = false;
      int i = 1;
      while (!found) {
         UTIL_CHECK(i <= nSpecies_);
         if (firstMoleculeId_[i] > moleculeId) {
            found = true;
         } else {
            ++i;
         }
      }
      --i;
      int diff = moleculeId - firstMoleculeId_[i];
      UTIL_CHECK(diff >= 0);
      int nMolecule = nMolecule_[i];
      UTIL_CHECK(diff < nMolecule);
      context.speciesId  = i;
      context.moleculeId = diff;
   }

   void SpeciesFinder::findPart(int partId, 
                                SpeciesFinder::Part& context)
   {
      UTIL_CHECK(nSpecies_ > 0);
      UTIL_CHECK(partId >= 0);
      UTIL_CHECK(partId < firstPartId_[nSpecies_]);
      bool found = false;
      int i = 1;
      while (!found) {
         UTIL_CHECK(i <= nSpecies_);
         if (firstPartId_[i] > partId) {
            found = true;
         } else {
            ++i;
         }
      }
      --i;
      int diff = partId - firstPartId_[i];
      UTIL_CHECK(diff >= 0);
      int nMolecule = nMolecule_[i];
      UTIL_CHECK(diff < nMolecule*nPart_[i]);
      context.speciesId  = i;
      context.moleculeId  = diff/nMolecule;
      context.partId  = diff - nMolecule*context.moleculeId;
   }

} 
