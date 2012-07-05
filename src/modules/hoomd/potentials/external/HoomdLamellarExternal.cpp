#ifdef INTER_EXTERNAL
#ifndef HOOMD_LAMELLAR_EXTERNAL_CPP
#define HOOMD_LAMELLAR_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdLamellarExternal.h"


namespace McMd
{
   char classNameHoomdLamellar[] = "HoomdLamellarExternal";

   /**
   * Default constructor.
   */
   HoomdLamellarExternal::HoomdLamellarExternal()   
    : HoomdExternal< EvaluatorExternalPeriodic, gpu_compute_periodic_forces, classNameHoomdLamellar >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdLamellarExternal::HoomdLamellarExternal(const HoomdLamellarExternal& other)
    : HoomdExternal< EvaluatorExternalPeriodic, gpu_compute_periodic_forces, 
          classNameHoomdLamellar >(other)
   {
      externalParameter_ = other.externalParameter_;
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
      read<double>(in, "externalParameter", externalParameter_);
      read<double>(in, "interfaceWidth", width_);
      read<int>(in, "periodicity", periodicity_);

      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].x = __int_as_scalar(perpDirection_);
         params_[i].y = Scalar(prefactor_[i]*externalParameter_);
         params_[i].z = Scalar(width_);
         params_[i].w = __int_as_scalar(periodicity_);
      } 
   }

   /*
   * set external potential parameter
   */
   void HoomdLamellarExternal::setExternalParameter(double externalParameter)
   {
      externalParameter_ = externalParameter;
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].y = prefactor_[i]*externalParameter_;
      }
   }

   /* 
   * Get external potential interaction strength.
   */
   double HoomdLamellarExternal::externalParameter() const
   {
      return externalParameter_;
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
