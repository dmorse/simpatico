/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PeriodicExternal.h"
#include <iostream>

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   PeriodicExternal::PeriodicExternal() 
    : externalParameter_(),
      nWaveVectors_(),
      C_(),
      periodicity_(),
      interfaceWidth_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   {  setClassName("PeriodicExternal"); }
   
   /* 
   * Copy constructor.
   */
   PeriodicExternal::PeriodicExternal(const PeriodicExternal& other)
    : externalParameter_(other.externalParameter_),
      nWaveVectors_(other.nWaveVectors_),
      C_(other.C_),
      periodicity_(other.periodicity_),
      interfaceWidth_(other.interfaceWidth_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      prefactor_.allocate(nAtomType_);
      for (int i=0; i < nAtomType_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
      waveVectors_.allocate(nWaveVectors_);
      for (int j=0; j < Dimension; ++j) {
        for (int i=0; i < nWaveVectors_; ++i) {
          waveVectors_[i][j] = other.waveVectors_[i][j];
        }
      }
      phases_.allocate(nWaveVectors_);
      for (int i=0; i < nWaveVectors_; ++i) {
        phases_[i] = other.phases_[i];
      }
      for (int i=0; i < Dimension; ++i) {
        shift_[i] = other.shift_[i];
      }
   } 
     
   /* 
   * Assignment operator.
   */
   PeriodicExternal& PeriodicExternal::operator = (const PeriodicExternal& other)
   {
      externalParameter_   = other.externalParameter_;
      nWaveVectors_        = other.nWaveVectors_;
      C_                   = other.C_;
      periodicity_         = other.periodicity_;
      interfaceWidth_      = other.interfaceWidth_;
      boundaryPtr_         = other.boundaryPtr_;
      nAtomType_           = other.nAtomType_;
      isInitialized_       = other.isInitialized_;
      for (int i=0; i < nAtomType_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
      for (int j=0; j < Dimension; ++j) {
        for (int i=0; i < nWaveVectors_; ++i) {
          waveVectors_[i][j] = other.waveVectors_[i][j];
        }
      }
      phases_.allocate(nWaveVectors_);
      for (int i=0; i < nWaveVectors_; ++i) {
        phases_[i] = other.phases_[i];
      }
      for (int i=0; i < Dimension; ++i) {
        shift_[i] = other.shift_[i];
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void PeriodicExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > PeriodicExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void PeriodicExternal::setExternalParameter(double externalParameter) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("PeriodicExternal potential is not initialized");
      }
     
      externalParameter_ = externalParameter;
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void PeriodicExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void PeriodicExternal::readParameters(std::istream &in) 
   {
      UTIL_CHECK(nAtomType_ > 0);
      UTIL_CHECK(boundaryPtr_);
  
      // Read parameters
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);
      read<double>(in, "externalParameter", externalParameter_);
      read<int>(in, "nWaveVectors", nWaveVectors_);
      read<double>(in, "C", C_);
      waveVectors_.allocate(nWaveVectors_);
      readDArray<Vector>(in, "waveVectors", waveVectors_, nWaveVectors_);
      phases_.allocate(nWaveVectors_);
      readDArray<double>(in, "phases", phases_, nWaveVectors_);
      read<Vector>(in, "shift", shift_);
      read<double>(in, "interfaceWidth", interfaceWidth_);
      read<int>(in, "periodicity", periodicity_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void PeriodicExternal::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nAtomType_ > 0);
      prefactor_.allocate(nAtomType_);
      loadDArray<double>(ar, "prefactor", prefactor_, nAtomType_);
      loadParameter<double>(ar, "externalParameter", externalParameter_);
      loadParameter<int>(ar, "nWavevectors", nWaveVectors_);
      loadParameter<double>(ar, "C", C_);
      waveVectors_.allocate(nWaveVectors_);
      loadDArray<Vector>(ar, "waveVectors", waveVectors_, nWaveVectors_);
      phases_.allocate(nWaveVectors_);
      loadDArray<double>(ar, "phases", phases_, nWaveVectors_);
      loadParameter<Vector>(ar, "shift", shift_);
      loadParameter<double>(ar, "interfaceWidth", interfaceWidth_);
      loadParameter<int>(ar, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void PeriodicExternal::save(Serializable::OArchive &ar)
   {
      ar << prefactor_;
      ar << externalParameter_;
      ar << nWaveVectors_;
      ar << C_;
      ar << waveVectors_;
      ar << phases_;
      ar << shift_;
      ar << interfaceWidth_;
      ar << periodicity_;
   }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void PeriodicExternal::set(std::string name, double value)
   {
      if (name == "externalParameter") {
         externalParameter_ = value;
      } else
      if (name == "C") {
         C_ = value;
      } else
      if (name == "interfaceWidth") {
         interfaceWidth_ = value;
      } else
      if (name == "periodicity") {
         periodicity_ = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double PeriodicExternal::get(std::string name) const
   {
      double value = 0.0;
      if (name == "externalParameter") {
         value = externalParameter_;
      } else
      if (name == "C") {
         value = C_;
      } else
      if (name == "interfaceWidth") {
         value = interfaceWidth_;
      } else
      if (name == "periodicity") {
         value = periodicity_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   double PeriodicExternal::externalParameter() const
   { 
     return externalParameter_; 
   }

   /*
   * Return name string "PeriodicExternal".
   */
   std::string PeriodicExternal::className() const
   {  return std::string("PeriodicExternal"); }
 
} 
