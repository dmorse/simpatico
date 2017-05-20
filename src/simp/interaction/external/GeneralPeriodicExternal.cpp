/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GeneralPeriodicExternal.h"

#include <iostream>

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   GeneralPeriodicExternal::GeneralPeriodicExternal() 
    : externalParameter_(),
      nWaveVectors_(),
      periodicity_(),
      interfaceWidth_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   { setClassName("GeneralPeriodicExternal"); }
   
   /* 
   * Copy constructor.
   */
   GeneralPeriodicExternal::GeneralPeriodicExternal(const GeneralPeriodicExternal& other)
    : externalParameter_(other.externalParameter_),
      nWaveVectors_(other.nWaveVectors_),
      periodicity_(other.periodicity_),
      interfaceWidth_(other.interfaceWidth_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      prefactor_.allocate(nAtomType_ * nWaveVectors_);
      for (int i=0; i < nAtomType_ * nWaveVectors_; ++i) {
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
      prefactor_.allocate(Dimension * nWaveVectors_);
      for (int i=0; i < Dimension * nWaveVectors_; ++i) {
        shifts_[i] = other.shifts_[i];
      }
   } 
     
   /* 
   * Assignment operator.
   */
   GeneralPeriodicExternal& GeneralPeriodicExternal::operator = (const GeneralPeriodicExternal& other)
   {
      externalParameter_   = other.externalParameter_;
      nWaveVectors_        = other.nWaveVectors_;
      periodicity_         = other.periodicity_;
      interfaceWidth_      = other.interfaceWidth_;
      boundaryPtr_         = other.boundaryPtr_;
      nAtomType_           = other.nAtomType_;
      isInitialized_       = other.isInitialized_;
      for (int i=0; i < nAtomType_ * nWaveVectors_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
      for (int j=0; j < Dimension; ++j) {
        for (int i=0; i < nWaveVectors_; ++i) {
          waveVectors_[i][j] = other.waveVectors_[i][j];
        }
      }
      for (int i=0; i < nWaveVectors_; ++i) {
        phases_[i] = other.phases_[i];
      }
      for (int i=0; i < Dimension * nWaveVectors_; ++i) {
        shifts_[i] = other.shifts_[i];
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void GeneralPeriodicExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > GeneralPeriodicExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void GeneralPeriodicExternal::setExternalParameter(double externalParameter) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("GeneralPeriodicExternal potential is not initialized");
      }
     
      externalParameter_ = externalParameter;
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void GeneralPeriodicExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void GeneralPeriodicExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
  
      // Read parameters
      read<int>(in, "nWaveVectors", nWaveVectors_);

      prefactor_.allocate(nAtomType_ * nWaveVectors_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_ * nWaveVectors_);

      waveVectors_.allocate(nWaveVectors_);
      readDArray<Vector>(in, "waveVectors", waveVectors_, nWaveVectors_);

      phases_.allocate(nWaveVectors_);
      readDArray<double>(in, "phases", phases_, nWaveVectors_);

      shifts_.allocate(Dimension * nWaveVectors_);
      readDArray<double>(in, "shifts", shifts_, Dimension * nWaveVectors_);

      read<double>(in, "externalParameter", externalParameter_);

      read<double>(in, "interfaceWidth", interfaceWidth_);

      read<int>(in, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void GeneralPeriodicExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_;
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      prefactor_.allocate(nAtomType_ * nWaveVectors_);
      loadDArray<double>(ar, "prefactor", prefactor_, nAtomType_ * nWaveVectors_);

      waveVectors_.allocate(nWaveVectors_);
      loadDArray<Vector>(ar, "waveVectors", waveVectors_, nWaveVectors_);

      phases_.allocate(nWaveVectors_);
      loadDArray<double>(ar, "phases", phases_, nWaveVectors_);

      shifts_.allocate(Dimension * nWaveVectors_);
      loadDArray<double>(ar, "shifts", shifts_, Dimension * nWaveVectors_);

      loadParameter<double>(ar, "externalParameter", externalParameter_);

      loadParameter<double>(ar, "interfaceWidth", interfaceWidth_);

      loadParameter<int>(ar, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void GeneralPeriodicExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << prefactor_;
      ar << externalParameter_;
      ar << waveVectors_;
      ar << phases_;
      ar << shifts_;
      ar << interfaceWidth_;
      ar << periodicity_;
   }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void GeneralPeriodicExternal::set(std::string name, double value)
   {
      if (name == "externalParameter") {
         externalParameter_ = value;
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
   double GeneralPeriodicExternal::get(std::string name) const
   {
      double value = 0.0;
      if (name == "externalParameter") {
         value = externalParameter_;
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

   double GeneralPeriodicExternal::externalParameter() const
   { 
     return externalParameter_; 
   }

   /*
   * Return name string "GeneralPeriodicExternal".
   */
   std::string GeneralPeriodicExternal::className() const
   {  return std::string("GeneralPeriodicExternal"); }
 
} 
