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
   
   const int BennettsMethod::TagDerivative[2] = {1, 2};
   const int BennettsMethod::TagFermi[2] = {8, 9};

   /* 
   * Constructor.
   */
   BennettsMethod::BennettsMethod(System& system)
    : SystemDiagnostic<System>(system),
      shift_(0.0),
      lowerShift_(0.0),
      shifts_(),
      communicatorPtr_(0),
      myId_(-1),
      nProcs_(0),
      lowerId_(-1),
      upperId_(-1),
      nParameter_(0),
      myParam_(),
      lowerParam_(),
      upperParam_(),
      nSamplePerBlock_(1),
      myAccumulator_(),
      upperAccumulator_(),
      myArg_(0.0),
      lowerArg_(0.0),
      myFermi_(0.0),
      lowerFermi_(0.0),
      upperFermi_(0.0),
      outputFile_(),
      isInitialized_(false)
   {  
      setClassName("BennettsMethod");
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
      nParameter_ = system.perturbation().getNParameters(); 
      myParam_.allocate(nParameter_);
      lowerParam_.allocate(nParameter_);
      upperParam_.allocate(nParameter_);
   }

   /*
   * Read parameters and initialize.
   */
   void BennettsMethod::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);

      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
 
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
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

      //If nSamplePerBlock != 0, open an output file for block averages.
      if (myAccumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from archive.
   */
   void BennettsMethod::loadParameters(Serializable::IArchive &ar)
   {
      Diagnostic::loadParameters(ar);

      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
 
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         shifts_.allocate(nProcs_);
         shifts_[nProcs_-1] = 0.0;
         loadDArray<double>(ar, "shifts", shifts_, nProcs_-1);
         shift_ = shifts_[myId_];
      } else {
         loadParameter<double>(ar, "shift", shift_);
      }
      #else
      loadParameter<double>(ar, "shift", shift_);
      #endif
     
      ar & myAccumulator_; 
      ar & upperAccumulator_; 

      // ar & myId_;
      // ar & nProcs_;
      // ar & lowerId_;
      // ar & upperId_;
      // ar & nParameter_;
      // ar & myParam_;
      
      ar & lowerParam_;
      ar & upperParam_;
      ar & myArg_;
      ar & lowerArg_;
      ar & myFermi_;
      ar & lowerFermi_;
      ar & upperFermi_;


      //If nSamplePerBlock != 0, open an output file for block averages.
      if (myAccumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to archive.
   */
   void BennettsMethod::save(Serializable::OArchive &ar)
   { ar & *this; }
   
   void BennettsMethod::setup()
   { 
     if (!isInitialized_) {
         UTIL_THROW("Object is not initialized");
      }

      communicatorPtr_->Bcast((void*)&shifts_[0], nProcs_, MPI::DOUBLE, 0);
      if ( myId_ != 0 ) {
         lowerShift_ = shifts_[lowerId_];
      } else {}

      myAccumulator_.clear();
      upperAccumulator_.clear();
   }

   void BennettsMethod::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;

      int myPort, upperPort;
      MPI::Request requestFermi[2];
      MPI::Status  status;

      // Exchange perturbation parameters and differences
      if (myId_ != 0 && myId_ != nProcs_ - 1) {

         myPort = myId_%2;
         upperPort = upperId_%2;

         for (int i = 0; i < nParameter_; ++i) {
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

         outputFile_ << Dbl(myFermi_) << "    " << Dbl(upperFermi_) << "    ";
         outputFile_ << std::endl;
 
         } else if (myId_ == 0) {
            
            myPort = myId_%2;
            upperPort = upperId_%2;

            for (int i = 0; i < nParameter_; ++i) {
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
            
            outputFile_ << Dbl(myFermi_) << "    " << Dbl(upperFermi_) << "    ";
            outputFile_ << std::endl;

         } else if (myId_ == nProcs_ - 1) {
            
            myPort = myId_%2;

            for (int i = 0; i < nParameter_; ++i) {
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
      double EnergyDiff, ratio, improvedShift;
      if (myId_ != nProcs_ - 1) {
       
         ratio = upperAccumulator_.average()/myAccumulator_.average();
      
         EnergyDiff = log(ratio)+shift_;
      
         improvedShift = EnergyDiff;
            
         shifts_[myId_] = improvedShift; 
         shift_ = improvedShift; 
      } else {}
   
   }
   
   
   // Output results to file after simulation is completed.
   void BennettsMethod::output() 
   {
      if (myAccumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      if (myId_ != nProcs_ - 1) {
         myAccumulator_.output(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << std::endl;
         upperAccumulator_.output(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << std::endl;
      }
      outputFile_.close();

      analyze(); 

      communicatorPtr_->Gather((const void *) &shift_, 1, MPI::DOUBLE, (void *) &shifts_[0], 1, MPI::DOUBLE, 0);
      int i;
      if (myId_ == 0) {
         fileMaster().openOutputFile(outputFileName("_all.dat"), outputFile_);
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
