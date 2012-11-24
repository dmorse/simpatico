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
   Perturbation::Perturbation(int size, int rank)
    : ParamComposite(),
      size_(size),
      rank_(rank),
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
      read<int>(in, "mode", mode_);
      read<int>(in, "nParameters", nParameters_);
      parameters_.allocate(size_, nParameters_);
      parameter_.allocate(nParameters_);
      initialParameter_.allocate(nParameters_);
      finalParameter_.allocate(nParameters_);
      if (mode_ == 0) {
        readDMatrix<double>(in, "parameters", parameters_, size_, nParameters_);
        for (int i = 0; i < nParameters_; ++i) {
           parameter_[i] = parameters_(rank_,i);
        }
      } else if (mode_ == 1) {
        readDArray<double>(in, "initialParameter", initialParameter_, nParameters_);
        readDArray<double>(in, "finalParameter", finalParameter_, nParameters_);
        int i, j;
        for (i  = 0; i < nParameters_; ++i) {
           parameters_(0,i) = initialParameter_[i];
           parameters_(size_-1,i) = finalParameter_[i];
           for (j= 1; j < size_-1; ++j) {
              parameters_(j, i) = parameters_(0, i) 
                                + j*((parameters_(size_-1, i) - parameters_(0, i))/(size_-1));
           }
        }
        for (i  = 0; i < nParameters_; ++i) {
           parameter_[i] = parameters_(rank_, i);
        }
      }

      setParameter();  // Modify parameter of associated System
   }

   /*
   * Load internal state from an archive.
   */
   void Perturbation::loadParameters(Serializable::IArchive &ar)
   {  
      loadParameter<int>(ar, "mode", mode_);
      loadParameter<int>(ar, "nParameters", nParameters_);
      parameter_.allocate(nParameters_);
      parameters_.allocate(size_, nParameters_);
      initialParameter_.allocate(nParameters_);
      finalParameter_.allocate(nParameters_);
      if (mode_ == 0) {
        loadDMatrix<double>(ar, "parameters", parameters_, size_, nParameters_);
      } else if (mode_ == 1) {
        loadDArray<double>(ar, "initialParameter", initialParameter_, nParameters_);
        loadDArray<double>(ar, "finalParameter", finalParameter_, nParameters_);
        ar & parameters_;
      }
      ar & parameter_;
      setParameter();  // Modify parameter of associated System
   }

   /*
   * Save internal state to an archive.
   */
   void Perturbation::save(Serializable::OArchive &ar)
   {
      ar & mode_;
      ar & nParameters_;
      if (mode_ == 0) {
        ar & parameters_;
      } else if (mode_ == 1) {
        ar & initialParameter_;
        ar & finalParameter_;
        ar & parameters_;
      }
      ar & parameter_;
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

   /*
   * Get parameter i of system id.
   */
   double Perturbation::parameter(int i, int id)
   { return parameters_(id,i); }

   /*
   * Get number of parameters per System. 
   */
   int Perturbation::getNParameters() const
   { return nParameters_; }

}

#endif 
#endif  // ifdef  MCMD_PERTURB
