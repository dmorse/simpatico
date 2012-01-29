#ifdef MCMD_EXTERNAL
#ifndef HOOMD_LAMELLAR_EXTERNAL_CPP
#define HOOMD_LAMELLAR_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdLamellarExternal.h"


namespace McMd
{
   /**
   * Default constructor.
   */
   HoomdLamellarExternal::HoomdLamellarExternal()   
    : HoomdExternal< EvaluatorExternalLamellar, gpu_compute_lamellar_forces >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdLamellarExternal::HoomdLamellarExternal(const HoomdLamellarExternal& other)
    : HoomdExternal< EvaluatorExternalLamellar, gpu_compute_lamellar_forces >(other)
   {
      orderParameter_ = other.orderParameter_;
   }

   /**
   * read parameters from file
   */
   void HoomdLamellarExternal::readParam(std::istream &in)
   {
      // Read parameters
      read<int>(in, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);
      read<double>(in, "orderParameter", orderParameter_);
      read<double>(in, "interfaceWidth", width_);
      read<int>(in, "periodicity", periodicity_);

      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].x = __int_as_scalar(perpDirection_);
         params_[i].y = Scalar(prefactor_[i]*orderParameter_);
         params_[i].z = Scalar(width_);
         params_[i].w = __int_as_scalar(periodicity_);
      } 
      isInitialized_ = true;
   }

   /*
   * set external potential parameter
   */
   void HoomdLamellarExternal::setExternalParameter(double orderParameter)
   {
      orderParameter_ = orderParameter;
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].y = prefactor_[i]*orderParameter_;
      }
   }

   /* 
   * Get external potential interaction strength.
   */
   double HoomdLamellarExternal::externalParameter() const
   {
      return orderParameter_;
   }

   /*
   * return the class name
   */
   std::string HoomdLamellarExternal::className() const
   {
      return "HoomdLamellarExternal";
   }

}

#endif
#endif
