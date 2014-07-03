#ifndef INTER_SIMPLE_PERIODIC_EXTERNAL_CPP
#define INTER_SIMPLE_PERIODIC_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SimplePeriodicExternal.h"

#include <iostream>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   SimplePeriodicExternal::SimplePeriodicExternal() 
    : externalParameter_(),
      nWaveVectors_(),
      periodicity_(),
      interfaceWidth_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   { setClassName("SimplePeriodicExternal"); }
   
   /* 
   * Copy constructor.
   */
   SimplePeriodicExternal::SimplePeriodicExternal(const SimplePeriodicExternal& other)
    : externalParameter_(other.externalParameter_),
      nWaveVectors_(other.nWaveVectors_),
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
      waveIntVectors_.allocate(nWaveVectors_);
      for (int j=0; j < Dimension; ++j) {
        for (int i=0; i < nWaveVectors_; ++i) {
          waveIntVectors_[i][j] = other.waveIntVectors_[i][j];
        }
      }
   } 
     
   /* 
   * Assignment operator.
   */
   SimplePeriodicExternal& SimplePeriodicExternal::operator = (const SimplePeriodicExternal& other)
   {
      externalParameter_   = other.externalParameter_;
      nWaveVectors_        = other.nWaveVectors_;
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
          waveIntVectors_[i][j] = other.waveIntVectors_[i][j];
        }
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void SimplePeriodicExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > SimplePeriodicExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void SimplePeriodicExternal::setExternalParameter(double externalParameter) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("SimplePeriodicExternal potential is not initialized");
      }
     
      externalParameter_ = externalParameter;
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void SimplePeriodicExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void SimplePeriodicExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
  
      // Read parameters
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);

      read<double>(in, "externalParameter", externalParameter_);

      read<int>(in, "nWaveVectors", nWaveVectors_);
      waveIntVectors_.allocate(nWaveVectors_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWaveVectors_);
      
      read<double>(in, "interfaceWidth", interfaceWidth_);
      read<int>(in, "periodicity", periodicity_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void SimplePeriodicExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_;
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      prefactor_.allocate(nAtomType_);
      loadDArray<double>(ar, "prefactor", prefactor_, nAtomType_);
      loadParameter<double>(ar, "externalParameter", externalParameter_);
      waveIntVectors_.allocate(nWaveVectors_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWaveVectors_);
      loadParameter<double>(ar, "interfaceWidth", interfaceWidth_);
      loadParameter<int>(ar, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void SimplePeriodicExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << prefactor_;
      ar << externalParameter_;
      ar << waveIntVectors_;
      ar << interfaceWidth_;
      ar << periodicity_;
   }


   double SimplePeriodicExternal::externalParameter() const
   { 
     return externalParameter_; 
   }

   /*
   * Return name string "SimplePeriodicExternal".
   */
   std::string SimplePeriodicExternal::className() const
   {  return std::string("SimplePeriodicExternal"); }
 
} 
#endif
