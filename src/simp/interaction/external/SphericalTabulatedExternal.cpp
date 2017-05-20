/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SphericalTabulatedExternal.h"

#include <iostream>

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   SphericalTabulatedExternal::SphericalTabulatedExternal() 
    : rMax_(0.0),
      rMaxSq_(0.0),
      dr_(0.0),
      boundaryPtr_(0),
      nr_(0),
      nAtomType_(0),
      isInitialized_(false)
   {  setClassName("SphericalTabulatedExternal"); }
   
   /* 
   * Copy constructor.
   */
   SphericalTabulatedExternal::SphericalTabulatedExternal(const SphericalTabulatedExternal& 
                                                          other)
    : rMax_(other.rMax_),
      rMaxSq_(other.rMaxSq_),
      dr_(0.0),
      boundaryPtr_(other.boundaryPtr_),
      nr_(0),
      nAtomType_(other.nAtomType_),
      isInitialized_(false)
   {
      if (other.potentials_.isAllocated()) {
         UTIL_CHECK(other.potentials_.capacity() == nAtomType_);
         potentials_.allocate(nAtomType_);
         int i, j;
         for (i = 0; i < nAtomType_; ++i) {
            UTIL_CHECK(other.potentials_[i].isAllocated());
            UTIL_CHECK(other.potentials_[i].capacity() == nr_);
            potentials_[i].allocate(nr_);
            for (j = 0; j < nr_; ++j) {
               potentials_[i][j] = other.potentials_[i][j];
            }
         }
      }
      isInitialized_ = other.isInitialized_;
   }

   /* 
   * Assignment operator.
   */
   SphericalTabulatedExternal& 
   SphericalTabulatedExternal::operator = (const SphericalTabulatedExternal& other)
   {
      rMax_          = other.rMax_;
      rMaxSq_        = other.rMaxSq_;
      dr_            = other.dr_;
      boundaryPtr_   = other.boundaryPtr_;
      boundaryPtr_   = other.boundaryPtr_;
      nr_            = other.nr_;
      nAtomType_     = other.nAtomType_;
      if (other.potentials_.isAllocated()) {
         UTIL_CHECK(other.potentials_.capacity() == nAtomType_);
         if (!potentials_.isAllocated()) {
            potentials_.allocate(nAtomType_);
         }
         UTIL_CHECK(potentials_.capacity() == nAtomType_);
         int i, j;
         for (i = 0; i < nAtomType_; ++i) {
            UTIL_CHECK(other.potentials_[i].capacity() == nr_);
            if (!potentials_[i].isAllocated()) {
               potentials_[i].allocate(nr_);
            }
            UTIL_CHECK(potentials_[i].capacity() == nr_);
            for (j = 0; j < nr_; ++j) {
               potentials_[i][j] = other.potentials_[i][j];
            }
         }
      }
      isInitialized_ = other.isInitialized_;

      return *this;
   }

   /* 
   * Set nAtomType
   */
   void SphericalTabulatedExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > SphericalTabulatedExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   /* 
   * Set pointer to the Boundary.
   */
   void SphericalTabulatedExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }

   /* 
   * Read potential parameters from file.
   */
   void SphericalTabulatedExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
   
      // Read parameters
      read<double>(in, "rMax", rMax_);
      read<int>(in, "nr", nr_);
      read<std::string>(in, "filename", filename_);

      rMaxSq_ = rMax_*rMax_;
      dr_ = rMax_/double(nr_ - 1);

      // Allocate arrays containing tabulated potentials
      potentials_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         potentials_[i].allocate(nr_);
      }
      std::cout << "Allocated potentials" << std::endl;

      // Read tabulated potential from file
      std::ifstream file;
      file.open(filename_.c_str());
      UTIL_CHECK(file.is_open());
      int i, j, k;
      for (j = 0; j < nr_; ++j) {
         file >> k;
         UTIL_CHECK(k == j);
         for (i = 0; i < nAtomType_; ++i) {
            file >> potentials_[i][j];
         }
      }
      file.close();

      // Shift potential to be zero at rMax
      std::cout << std::endl;
      for (j = 0; j < nr_; ++j) {
         std::cout << "  " << j;
         for (i = 0; i < nAtomType_; ++i) {
            potentials_[i][j] -= potentials_[i][nr_ - 1];
            std::cout << "  " << potentials_[i][j];
         }
         std::cout << std::endl;
      }
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void SphericalTabulatedExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_; 
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before loadParameters");
      }
      loadParameter<double> (ar, "rMax", rMax_);
      loadParameter<int> (ar, "nr", nr_);
      rMaxSq_ = rMax_*rMax_;
      dr_ = rMax_/double(nr_ - 1);

      // Allocate and read arrays containing tabulated potentials
      if (!potentials_.isAllocated()) {
         potentials_.allocate(nAtomType_);
      }
      UTIL_CHECK(potentials_.capacity() == nAtomType_);
      int i, j;
      for (i = 0; i < nAtomType_; ++i) {
         if (!potentials_[i].isAllocated()) {
            potentials_[i].allocate(nr_);
         }
         UTIL_CHECK(potentials_[i].capacity() == nr_);
         for (j = 0; j < nr_; ++j) {
            ar >> potentials_[i][j];
         }
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void SphericalTabulatedExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << rMax_;
      ar << nr_;
      UTIL_CHECK(potentials_.isAllocated());
      int i, j;
      for (i = 0; i < nAtomType_; ++i) {
         UTIL_CHECK(potentials_[i].isAllocated());
         for (j = 0; j < nr_; ++j) {
            ar << potentials_[i][j];
         }
      }
   }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void SphericalTabulatedExternal::set(std::string name, double value)
   {
      if (name == "rMax") {
         rMax_ = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double SphericalTabulatedExternal::get(std::string name) const
   {
      double value = 0.0;
      if (name == "rMax") {
         value = rMax_;
      } else 
      if (name == "nr") {
         value = nr_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /**
   * Return name string "SphericalTabulatedExternal" for this interaction class.
   */
   std::string SphericalTabulatedExternal::className() const
   {  return std::string("SphericalTabulatedExternal"); }
 

} 
