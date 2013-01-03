#ifndef INTER_MULTI_HARMONIC_DIHEDRAL_CPP
#define INTER_MULTI_HARMONIC_DIHEDRAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MultiHarmonicDihedral.h"

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MultiHarmonicDihedral::MultiHarmonicDihedral()
    : nDihedralType_(0)
   { 
      setClassName("MultiHarmonicDihedral");
   }

   /* 
   * Copy constructor.
   */
   MultiHarmonicDihedral::MultiHarmonicDihedral(const MultiHarmonicDihedral& other)
    : nDihedralType_(other.nDihedralType_)
   { for (int i = 0; i < nDihedralType_; ++i) kappa_[i] = other.kappa_[i]; }

   /* 
   * Assignment.
   */
   MultiHarmonicDihedral& MultiHarmonicDihedral::operator = (const MultiHarmonicDihedral& other)
   {
      nDihedralType_ = other.nDihedralType_;
      for (int i = 0; i < nDihedralType_; ++i) kappa_[i] = other.kappa_[i];
      return *this;
   }

   /* 
   * Set the nDihedralType_ member.
   */
   void MultiHarmonicDihedral::setNDihedralType(int nDihedralType)
   {  
      if (nDihedralType <= 0) {
         UTIL_THROW("nDihedralType must be positive");
      }
      nDihedralType_ = nDihedralType;
      parameters.allocate(nDihedralType_);
   }

   /* 
   * Read bend interaction parameters kappa from file.
   */
   void MultiHarmonicDihedral::readParameters(std::istream &in) 
   {
      if (nDihedralType_ <= 0) {
         UTIL_THROW("nDihedralType must be set before readParam");
      }
      read<ParameterSet>(in, "K",  parameters[i]);
   }

   /*
   * Load internal state from an archive.
   */
   void MultiHarmonicDihedral::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nDihedralType_; 
      if (nDihedralType_ <= 0) {
         UTIL_THROW( "nDihedralType must be positive");
      }
      parameters.allocate(nDihedralType_);
      for (int i = 0; i < nDihedralType_; ++i) {
         ParametersSet* p = &parameters[i];
         ar & p->K0;
         ar & p->K1;
         ar & p->K2;
         ar & p->K3;
         ar & p->K4;
         ar & p->A0;
         ar & p->A1;
         ar & p->A2;
         ar & p->A3;
         ar & p->A4;
      }
   }

   /*
   * Save internal state to an archive.
   */
   void MultiHarmonicDihedral::save(Serializable::OArchive &ar)
   {
      ar << nDihedralType_;
      for (int i = 0; i < nDihedralType_; ++i) {
         ParametersSet* p = &parameters[i];
         ar & p->K0;
         ar & p->K1;
         ar & p->K2;
         ar & p->K3;
         ar & p->K4;
         ar & p->A0;
         ar & p->A1;
         ar & p->A2;
         ar & p->A3;
         ar & p->A4;
      }
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void MultiHarmonicDihedral::set(std::string name, int type, double value)
   {
      if (name == "K0") {
        parameters[type].K0 = value;
      } else 
      if (name == "K1") {
        parameters[type].K1 = value;
      } else 
      if (name == "K2") {
        parameters[type].K2 = value;
      } else 
      if (name == "K3") {
        parameters[type].K3 = value;
      } else 
      if (name == "K4") {
        parameters[type].K4 = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double MultiHarmonicDihedral::get(std::string name, int type) const
   {
      double value = 0.0;
      ParameterSet p = &parameters[type];
      if (name == "K0") {
         value = p->K0;
      } else 
      if (name == "K1") {
         value = p->K1;
      } else 
      if (name == "K2") {
         value = p->K2;
      } else 
      if (name == "K3") {
         value = p->K3;
      } else 
      if (name == "K4") {
         value = p->K4;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /*
   * Return name string "MultiHarmonicDihedral" for this evaluator class.
   */
   std::string MultiHarmonicDihedral::className() const
   {  return std::string("MultiHarmonicDihedral"); }

} 
#endif
