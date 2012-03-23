#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_BENNETTS_METHOD_CPP
#define MCMD_BENNETTS_METHOD_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BennettsMethod.h"           // class header
#include <mcMd/perturb/Perturbation.h>  
#include <mcMd/perturb/LinearPerturbation.h>  
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/format/Dbl.h>
#ifdef UTIL_MPI
#include <util/containers/DArray.h>     // member
#endif

#include <cstdio> 
#include <sstream> 

namespace McMd
{

   using namespace Util;
   
   const int BennettsMethod::TagDerivative[2]    = {1, 2};
   const int BennettsMethod::TagFermi[2]    = {8, 9};

   /* 
   * Constructor.
   */
   BennettsMethod::BennettsMethod(System& system)
    : SystemDiagnostic<System>(system),
      communicatorPtr_(0),
      myId_(-1),
      nProcs_(0),
      lowerId_(-1),
      upperId_(-1),
      nParameters_(0),
      shift_(0.0),
      myAccumulator_(),
      upperAccumulator_(),
      nSamplePerBlock_(1),
      lowerShift_(0.0),
      myParam_(),
      lowerParam_(),
      upperParam_(),
      myDerivative_(0.0),
      myArg_(0.0),
      lowerArg_(0.0),
      myFermi_(0.0),
      lowerFermi_(0.0),
      upperFermi_(0.0),
      outputFile_()
   {  
      communicatorPtr_ = &(system.simulation().communicator());
      myId_   = communicatorPtr_->Get_rank();
      nProcs_ = communicatorPtr_->Get_size();   
      if (myId_ != 0) { 
         lowerId_ = myId_ - 1;
      } else {
         lowerId_ = myId_;
      }
      if (myId_ != nProcs_-1) {
         upperId_ = myId_ + 1;
      } else {
         upperId_ = myId_;
      }
      nParameters_ = system.perturbation().getNParameters(); 
      myParam_.allocate(nParameters_);
      lowerParam_.allocate(nParameters_);
      upperParam_.allocate(nParameters_);
   }

   /*
   * Read parameters and initialize.
   */
   void BennettsMethod::readParam(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);

      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
 
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         shifts_.allocate(nProcs_);
         shifts_[nProcs_-1]=0.0;
         readDArray<double>(in, "shifts", shifts_, nProcs_-1);
         shift_ = shifts_[myId_];
      } else {
         read<double>(in, "shift", shift_);
      }
      #else
      read<double>(in, "shift", shift_);
      #endif
      
      myAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      upperAccumulator_.setNSamplePerBlock(nSamplePerBlock_);

   }
    
   void BennettsMethod::initialize()
   { 
      communicatorPtr_->Bcast((void*)&shifts_[0], nProcs_, MPI::DOUBLE, 0);
      if ( myId_ != 0 ) {
         lowerShift_ = shifts_[lowerId_];
      } else {}
   }

   void BennettsMethod::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;

      int myPort, upperPort;
      MPI::Request requestFermi[2];
      MPI::Status  status;
      // Exchange derivatives with partner
      if (myId_ != 0 && myId_ != nProcs_ - 1) {
         myPort = myId_%2;
         upperPort = upperId_%2;

         for (int i = 0; i < nParameters_; ++i) {
            myParam_[i] = system().perturbation().parameter(i);
            lowerParam_[i] = system().perturbation().parameter(i,lowerId_);
            upperParam_[i] = system().perturbation().parameter(i,upperId_);
         }
         
         myArg_ = system().perturbation().difference(upperParam_);
         lowerArg_ = system().perturbation().difference(lowerParam_);

         myArg_ -= shift_;
         lowerArg_ += lowerShift_;
            
         myFermi_ = 1/(1 + exp(myArg_));
         lowerFermi_ = 1/(1 + exp(lowerArg_));
            
         requestFermi[0] = communicatorPtr_->Irecv(&upperFermi_, 1, MPI::DOUBLE, upperId_,
                                             TagFermi[upperPort]);
            
         requestFermi[1] = communicatorPtr_->Isend(&lowerFermi_, 1, MPI::DOUBLE, lowerId_,
                                             TagFermi[myPort]);
            
         // Synchronizing
         requestFermi[0].Wait();
         requestFermi[1].Wait();

         upperAccumulator_.sample(upperFermi_);
         myAccumulator_.sample(myFermi_);
 
         } else if (myId_ == 0) {
            
            myPort = myId_%2;
            upperPort = upperId_%2;

            for (int i = 0; i < nParameters_; ++i) {
               myParam_[i] = system().perturbation().parameter(i);
               upperParam_[i] = system().perturbation().parameter(i, upperId_);
            }

            myArg_ = system().perturbation().difference(upperParam_);
            myArg_ -= shift_;
            
            myFermi_ = 1/(1 + exp(myArg_));
            
            requestFermi[0] = communicatorPtr_->Irecv(&upperFermi_, 1, MPI::DOUBLE, upperId_,
                                              TagFermi[upperPort]);
             
            // Synchronizing
            requestFermi[0].Wait();
            
            upperAccumulator_.sample(upperFermi_);
            myAccumulator_.sample(myFermi_);

         } else if (myId_ == nProcs_ - 1) {
            
            myPort = myId_%2;

            for (int i = 0; i < nParameters_; ++i) {
               myParam_[i] = system().perturbation().parameter(i);
               lowerParam_[i] = system().perturbation().parameter(i, lowerId_);
            }

            lowerArg_ = system().perturbation().difference(lowerParam_);
            lowerArg_ += lowerShift_;
            
            lowerFermi_ = 1/(1 + exp(lowerArg_));
            
            requestFermi[1] = communicatorPtr_->Isend(&lowerFermi_, 1, MPI::DOUBLE, lowerId_,
                                               TagFermi[myPort]);
         
            // Synchronizing
            requestFermi[1].Wait();
        }
   }

   void BennettsMethod::analyze()
   {
      double EnergyDiff_, ratio_, improvedShift_;
      if (myId_ != nProcs_ - 1) {
       
         ratio_ = upperAccumulator_.average()/myAccumulator_.average();
      
         EnergyDiff_ = log(ratio_)+shift_;
      
         improvedShift_ = EnergyDiff_;
            
         shifts_[myId_] = improvedShift_; 
         shift_ = improvedShift_; 
    
      } else {}
   
   }
   
   
   // Output results to file after simulation is completed.
   
   void BennettsMethod::output() 
   {
      analyze(); 

      communicatorPtr_->Gather((const void *) &shift_, 1, MPI::DOUBLE, (void *) &shifts_[0], 1, MPI::DOUBLE, 0);
      int i;
      if (myId_ == 0) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
         for (i = 0; i < nProcs_-1; ++i) {
            outputFile_ << Dbl(shifts_[i]);
            outputFile_ << std::endl;
         }
         outputFile_.close();
      
      } else {}
   }
}
#endif    // ifndef BENNETTS_METHOD_CPP
#endif    // ifdef  UTIL_MPI
#endif    // ifdef  MCMD_PERTURB
