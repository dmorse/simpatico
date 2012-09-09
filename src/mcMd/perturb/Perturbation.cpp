#ifdef  MCMD_PERTURB
#ifndef MCMD_PERTURBATION_CPP
#define MCMD_PERTURBATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Perturbation.h"  // class header

namespace McMd
{

   /*
   * Constructor.
   */
   Perturbation::Perturbation()
    : ParamComposite(),
      nParameters_(0),
      mode_(0)
   {  setClassName("Perturbation"); }

   /*
   * Destructor.
   */
   Perturbation::~Perturbation()
   {}

   /*
   * Read parameter(s) and set system parameters.
   */
   void Perturbation::readParameters(std::istream& in)
   {  
      #ifdef UTIL_MPI
      int i, j;
      if (hasParamCommunicator()) {
         int size = paramCommunicator().Get_size();
         int rank = paramCommunicator().Get_rank();
         read<int>(in, "mode", mode_);
         read<int>(in, "nParameters", nParameters_);
         parameters_.allocate(size, nParameters_);
         parameter_.allocate(nParameters_);
         initialParameter_.allocate(nParameters_);
         finalParameter_.allocate(nParameters_);
         if (mode_ == 0) {
           readDMatrix<double>(in, "parameters", parameters_, size, nParameters_);
           for (i = 0; i < nParameters_; ++i) {
              parameter_[i] = parameters_(rank,i);
           }
         } else if (mode_ == 1) {
           readDArray<double>(in, "initialParameter", initialParameter_, nParameters_);
           readDArray<double>(in, "finalParameter", finalParameter_, nParameters_);
           for (i  = 0; i < nParameters_; ++i) {
              parameters_(0,i) = initialParameter_[i];
              parameters_(size-1,i) = finalParameter_[i];
              for (j= 1; j < size-1; ++j) {
                 parameters_(j, i) = parameters_(0, i) 
                                   + j*((parameters_(size-1, i) - parameters_(0, i))/(size-1));
              }
           }
           for (i  = 0; i < nParameters_; ++i) {
              parameter_[i] = parameters_(rank, i);
           }
         }
      } else {
         parameter_.allocate(nParameters_);
         readDArray<double>(in, "parameter", parameter_, nParameters_);
      }
      #else
      parameter_.allocate(nParameters_);
      readDArray<double>(in, "parameter", parameter_, nParameters_);
      #endif

      setParameter();  // Modify parameter of associated System
   }

   /*
   * Set the perturbation parameter.
   *
   * Set the parameter of the associated system.
   */
   void Perturbation::setParameter(DArray<double> parameter)
   {
      parameter_ = parameter;   // Set the class member.
      setParameter();           // Modify associated System.
   }
   
   #ifdef UTIL_MPI
   double Perturbation::parameter(int i, int id)
   { return parameters_(id,i); }
   #endif

   int Perturbation::getNParameters() const
   { return nParameters_; }

}

#endif 
#endif  // ifdef  MCMD_PERTURB
