#ifdef SIMP_EXTERNAL
#ifndef HOOMD_PERIODIC_EXTERNAL_CPP
#define HOOMD_PERIODIC_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdPeriodicExternal.h"


namespace McMd
{
   char classNameHoomdPeriodic[] = "HoomdPeriodicExternal";

   /**
   * Default constructor.
   */
   HoomdPeriodicExternal::HoomdPeriodicExternal()   
    : HoomdExternal< EvaluatorPeriodicExternal, gpu_compute_periodic_forces, classNameHoomdPeriodic >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdPeriodicExternal::HoomdPeriodicExternal(const HoomdPeriodicExternal& other)
    : HoomdExternal< EvaluatorPeriodicExternal, gpu_compute_periodic_forces, 
          classNameHoomdPeriodic >(other)
   {
      externalParameter_ = other.externalParameter_;
      interfaceWidth_ = other.interfaceWidth_;
      periodicity_ = other.periodicity_;
   }

   /**
   * read parameters from file
   */
   void HoomdPeriodicExternal::readParameters(std::istream &in)
   {
      // Read parameters
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);

      read<double>(in, "externalParameter", externalParameter_);

      waveIntVectors_.allocate(3);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, 3);
      read<double>(in, "interfaceWidth", interfaceWidth_);
      read<int>(in, "periodicity", periodicity_);
   
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].order_parameter = Scalar(prefactor_[i]*externalParameter_);
         params_[i].lattice_vector_1 = make_int3(waveIntVectors_[0][0], waveIntVectors_[0][1], waveIntVectors_[0][2]);
         params_[i].lattice_vector_2 = make_int3(waveIntVectors_[1][0], waveIntVectors_[1][1], waveIntVectors_[2][2]);
         params_[i].lattice_vector_3 = make_int3(waveIntVectors_[2][0], waveIntVectors_[2][1], waveIntVectors_[2][2]);
         params_[i].interface_width = Scalar(interfaceWidth_);
         params_[i].periodicity = periodicity_;
      }
   }

   /*
   * set external potential parameter
   */
   void HoomdPeriodicExternal::setExternalParameter(double externalParameter)
   {
      externalParameter_ = externalParameter;
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].order_parameter = prefactor_[i]*externalParameter_;
      }
   }

   /* 
   * Get external potential interaction strength.
   */
   double HoomdPeriodicExternal::externalParameter() const
   {
      return externalParameter_;
   }

   /*
   * return the class name
   */
   std::string HoomdPeriodicExternal::className() const
   {
      return "HoomdPeriodicExternal";
   }

}

#endif
#endif
