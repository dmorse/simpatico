#ifndef BENNETS_METHOD_CPP
#define MCMD_BENNETS_METHOD_CPP

#ifdef UTIL_MPI
#ifdef MCMD_PERTURB

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BennetsMethod.h"           // class header
#include <mcMd/perturb/Perturbation.h>  
#include <mcMd/perturb/LinearPerturbation.h>  
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <util/mpi/MpiSendRecv.h>
#ifdef UTIL_MPI
#include <util/containers/DArray.h>     // member
#endif

#include <cstdio> 
#include <sstream> 

namespace McMd
{

   using namespace Util;
   
   const int BennetsMethod::TagDerivative[2]    = {1, 2};
   const int BennetsMethod::TagFermi[2]    = {8, 9};

   /* 
   * Constructor.
   */
   BennetsMethod::BennetsMethod(System& system)
    : SystemDiagnostic<System>(system),
      communicatorPtr_(0),
      nProcs_(0),
      myId_(-1),
      lowerId_(-1),
      upperId_(-1),
      shift_(0.0),
      lowerShift_(0.0),
      myParam_(0.0),
      lowerParam_(0.0),
      upperParam_(0.0),
      myDerivative_(0.0),
      upperDerivative_(0.0),
      myArg_(0.0),
      lowerArg_(0.0),
      myFermi_(0.0),
      upperFermi_(0.0),
      lowerFermi_(0.0),
      myOutputFile_(),
      upperOutputFile_(),
      myAccumulator_(),
      upperAccumulator_(),
      nSamplePerBlock_(1),
      logger()
   {  
      communicatorPtr_ = &(system.simulation().communicator());
      myId_   = communicatorPtr_->Get_rank();
      nProcs_ = communicatorPtr_->Get_size();
      myParam_ = system.perturbation().parameter(myId_);
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
      lowerParam_ = system.perturbation().parameter(lowerId_); 
      upperParam_ = system.perturbation().parameter(upperId_); 
   }

   /*
   * Read parameters and initialize.
   */
   void BennetsMethod::readParam(std::istream& in)
   {
      readInterval(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
 
      #ifdef UTIL_MPI
      if (hasParamCommunicator()) {
         shifts_.allocate(nProcs_);
         readDArray<double>(in, "shifts", shifts_, nProcs_);
         shift_ = shifts_[myId_];
         
      } else {
         read<double>(in, "shift", shift_);
      }
      #else
      read<double>(in, "shift", shift_);
      #endif
      

      myAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      upperAccumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (myId_ != nProcs_ -1) {
         if (myAccumulator_.nSamplePerBlock()) {
            
            fileMaster().openOutputFile("current.dat", myOutputFile_);
         
         }
         if (upperAccumulator_.nSamplePerBlock()) {
             
            fileMaster().openOutputFile("upper.dat", upperOutputFile_);
      
        }
      } else {}
   }
    
 
   void BennetsMethod::Shifts()
   { 
     int myPort, upperPort; 
     MPI::Request requestDerivative[2];
     MPI::Status  status;
     
     if (myId_ != 0 && myId_ != nProcs_ - 1) {

         myPort =  myId_%2;
         upperPort = upperId_%2;

         myDerivative_ = system().perturbation().derivative();
         
         requestDerivative[0] = communicatorPtr_->Irecv(&upperDerivative_, 1, MPI::DOUBLE, upperId_,
                                              TagDerivative[upperPort]);
         requestDerivative[1] = communicatorPtr_->Isend(&myDerivative_, 1, MPI::DOUBLE, lowerId_,
                                               TagDerivative[myPort]);

         // Synchronizing
         requestDerivative[0].Wait();
         requestDerivative[1].Wait();
         
         shifts_[myId_] = (0.5*(myDerivative_ + upperDerivative_)*(upperParam_ - myParam_)); 
         shift_ = shifts_[myId_];

     } else if (myId_ == 0) {

         upperPort = upperId_%2;
         
         myDerivative_ = system().perturbation().derivative();

         requestDerivative[0] = communicatorPtr_->Irecv(&upperDerivative_, 1, MPI::DOUBLE, upperId_,
                                              TagDerivative[upperPort]);

         // Synchronizing
         requestDerivative[0].Wait();

         shifts_[myId_] = (0.5*(myDerivative_ + upperDerivative_)*(upperParam_ - myParam_)); 
         shift_ = shifts_[myId_];

     } else if (myId_ == nProcs_ - 1) {

         myPort =  myId_%2;
         
         myDerivative_ = system().perturbation().derivative();

         requestDerivative[1] = communicatorPtr_->Isend(&myDerivative_, 1, MPI::DOUBLE, lowerId_,
                                               TagDerivative[myPort]);

         // Synchronizing
         requestDerivative[1].Wait();
         
     }
    
   }
 
   void BennetsMethod::initialize()
   { 
     
     Shifts();
     for (int i=0; i < (nProcs_ - 1); i++) {
           communicatorPtr_->Bcast(&shifts_[i], 1, MPI::DOUBLE, i);
     }
     if ( myId_ != 0 ) {
        lowerShift_ = shifts_[lowerId_];
     } else {}
   }

   void BennetsMethod::sample(long iStep) 
   {
         int myPort, upperPort;
         MPI::Request requestFermi[2];
         MPI::Status  status;
         // Exchange derivatives with partner
         if (myId_ != 0 && myId_ != nProcs_ - 1) {
            
            myPort = myId_%2;
            upperPort = upperId_%2;
            
            myDerivative_ = system().perturbation().derivative();

            myArg_ = (myDerivative_*(upperParam_ - myParam_))-shift_;
            lowerArg_ = (myDerivative_*(lowerParam_ - myParam_))+lowerShift_;
            
            myFermi_ = 1/(1 + exp(myArg_));
            lowerFermi_ = 1/(1 + exp(lowerArg_));
            
            requestFermi[0] = communicatorPtr_->Irecv(&upperFermi_, 1, MPI::DOUBLE, upperId_,
                                              TagFermi[upperPort]);
            
            requestFermi[1] = communicatorPtr_->Isend(&lowerFermi_, 1, MPI::DOUBLE, lowerId_,
                                               TagFermi[myPort]);
            
            // Synchronizing
            requestFermi[0].Wait();
            requestFermi[1].Wait();

            upperAccumulator_.sample(upperFermi_, upperOutputFile_);
            myAccumulator_.sample(myFermi_, myOutputFile_);
 
         } else if (myId_ == 0) {
            
            myPort = myId_%2;
            upperPort = upperId_%2;
            
            myDerivative_ = system().perturbation().derivative();
            
            myArg_ = (myDerivative_*(upperParam_ - myParam_))-shift_;
            
            myFermi_ = 1/(1 + exp(myArg_));
            
            requestFermi[0] = communicatorPtr_->Irecv(&upperFermi_, 1, MPI::DOUBLE, upperId_,
                                              TagFermi[upperPort]);
             
            // Synchronizing
            requestFermi[0].Wait();
            
            upperAccumulator_.sample(upperFermi_, upperOutputFile_);
            myAccumulator_.sample(myFermi_, myOutputFile_);

         } else if (myId_ == nProcs_ - 1) {
            
            myPort = myId_%2;
            
            myDerivative_ = system().perturbation().derivative();
           
            lowerArg_ = (myDerivative_*(lowerParam_ - myParam_))+lowerShift_;
            
            lowerFermi_ = 1/(1 + exp(lowerArg_));
            
            requestFermi[1] = communicatorPtr_->Isend(&lowerFermi_, 1, MPI::DOUBLE, lowerId_,
                                               TagFermi[myPort]);
         
            // Synchronizing
            requestFermi[1].Wait();
        }
   }

   void BennetsMethod::analyze()
   {
     double EnergyDiff_, ratio_, improvedShift_;
     if (myId_ != nProcs_ - 1) {
       
        ratio_ = upperAccumulator_.average()/myAccumulator_.average();
      
        EnergyDiff_ = log(ratio_)+shift_;
      
        improvedShift_ = EnergyDiff_;
      
       //logger.begin(); 
       std::cout << "For processors " << myId_ << " and " << upperId_ << ", the free energy estimate is " << EnergyDiff_ << std::endl;
   
       std::cout << "For processors " << myId_ << " and " << upperId_ << ", the estimate for improved shift is " << improvedShift_ << std::endl;
    
       //logger.end();   
       shifts_[myId_] = improvedShift_; 
       shift_ = improvedShift_; 
    
    } else {}
   
   }
   
   
   // Output results to file after simulation is completed.
   
   void BennetsMethod::output() 
   {
      analyze(); 

      communicatorPtr_->Gather((const void *) &shift_, 1, MPI::DOUBLE, (void *) &shifts_[0], 1, MPI::DOUBLE, 0);

      /*if (myId_ == 1) {
+
         std::cout << "n10" << shifts_[0] << std::endl; 
         std::cout << "n10" << shifts_[1] << std::endl; 
         std::cout << "n10" << shifts_[2] << std::endl; 
      
      }*/
      
      // If outputFile_ was used to write block averages, close it.
      if (myId_ != nProcs_ -1) {
         if (myAccumulator_.nSamplePerBlock()) {
            
            myOutputFile_.close();
      
         }
      
         if (upperAccumulator_.nSamplePerBlock()) {
             
            upperOutputFile_.close();
      
         }
         
         fileMaster().openOutputFile("mysummary", myOutputFile_);
         writeParam(myOutputFile_); 
         myOutputFile_ << std::endl;
         myAccumulator_.output(myOutputFile_); 
         myOutputFile_.close();
      
         fileMaster().openOutputFile("uppersummary", upperOutputFile_);
         writeParam(upperOutputFile_); 
         upperOutputFile_ << std::endl;
         upperAccumulator_.output(upperOutputFile_); 
         upperOutputFile_.close();
      } else {}
   }
}
#endif    // ifdef  MCMD_PERTURB
#endif    // ifdef  UTIL_MPI
#endif    // ifndef BENNETS_METHOD_CPP
