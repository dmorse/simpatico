#ifdef SIMP_EXTERNAL
#ifndef HOOMD_LOCAL_LAMELLAR_ORDERING_EXTERNAL_CPP
#define HOOMD_LOCAL_LAMELLAR_ORDERING_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdLocalLamellarOrderingExternal.h"


namespace McMd
{
   char classNameHoomdLocalLamellarOrdering[] = "HoomdLocalLamellarOrderingExternal";

   /**
   * Default constructor.
   */
   HoomdLocalLamellarOrderingExternal::HoomdLocalLamellarOrderingExternal()   
    : HoomdExternal< EvaluatorLocalExternal, gpu_compute_local_forces, classNameHoomdLocalLamellarOrdering >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdLocalLamellarOrderingExternal::HoomdLocalLamellarOrderingExternal(const HoomdLocalLamellarOrderingExternal& other)
    : HoomdExternal< EvaluatorLocalExternal, gpu_compute_local_forces, 
          classNameHoomdLocalLamellarOrdering >(other)
   {
      externalParameter_ = other.externalParameter_;
   }

   /**
   * read parameters from file
   */
   void HoomdLocalLamellarOrderingExternal::readParameters(std::istream &in)
   {
      // Read parameters
      read<int>(in, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      read<int>(in, "parallelDirection", parallelDirection_);
      if (parallelDirection_ < 0 || parallelDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for parallel direction.");
      }
      read<double>(in, "fraction", fraction_);
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);
      read<double>(in, "externalParameter", externalParameter_);
      read<double>(in, "interfaceWidth", width_);
      read<int>(in, "periodicity", periodicity_);

      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].perp_index = perpDirection_;
         params_[i].parallel_index = parallelDirection_;
         params_[i].fraction_length = fraction_;
         params_[i].order_parameter = Scalar(prefactor_[i]*externalParameter_);
         params_[i].interface_width = Scalar(width_);
         params_[i].periodicity = periodicity_;
      } 
   }

   /*
   * set external potential parameter
   */
   void HoomdLocalLamellarOrderingExternal::setExternalParameter(double externalParameter)
   {
      externalParameter_ = externalParameter;
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].order_parameter = prefactor_[i]*externalParameter_;
      }
   }

   /* 
   * Get external potential interaction strength.
   */
   double HoomdLocalLamellarOrderingExternal::externalParameter() const
   {
      return externalParameter_;
   }

   /*
   * return the class name
   */
   std::string HoomdLocalLamellarOrderingExternal::className() const
   {
      return "HoomdLocalLamellarOrderingExternal";
   }

}

#endif
#endif
