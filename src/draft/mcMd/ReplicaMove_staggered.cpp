#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_REPLICA_MOVE_CPP
#define MCMD_REPLICA_MOVE_CPP


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ReplicaMove.h"                      

#include <mcMd/perturb/Perturbation.h>
#include <mcMd/perturb/LinearPerturbation.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <util/boundary/OrthorhombicBoundary.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/util/Observer.h>

#include <sstream>
#include <cmath>

namespace McMd
{

   using namespace Util;

   // Define static constant MPI message Tag.

   const int ReplicaMove::TagParam[2]    = {10, 11};

   const int ReplicaMove::TagDecision[2] = {20, 21};

   const int ReplicaMove::TagConfig[2]   = {30, 31};

   /*
   * Constructor.
   */
   ReplicaMove::ReplicaMove(System& system) :
      communicatorPtr_(0),
      myId_(-1),
      nProcs_(0),
      outputFile_(),
      repxAttempt_(),
      repxAccept_(),
      interval_(-1),
      nParameters_(0),
      myParam_(),
      ptParam_(),
      systemPtr_(&system),
      ptId_(-1),
      stepCount_(0),
      ptPositionPtr_(0),
      myPositionPtr_(0)
   {

      // Precondition
      if (!system.hasPerturbation()) {
         UTIL_THROW("Parent System has no Perturbation");
      }

      communicatorPtr_ = &(system.simulation().communicator());
      myId_   = communicatorPtr_->Get_rank();
      nProcs_ = communicatorPtr_->Get_size();

      // Generate output file name and open the file.
      //std::stringstream sMyId;
      //std::string fName("repx");
      //sMyId << myId_;
      //fName += sMyId.str();
      system.fileMaster().openOutputFile("repx", outputFile_);

      // Generate output file name and open the file.
      //std::stringstream energysMyId;
      //std::string energyfName("diff_energy");
      //energysMyId << myId_;
      //energyfName += energysMyId.str();

      // Initialize statistical accumulators.
      repxAttempt_[0] = 0;
      repxAttempt_[1] = 0;
      repxAccept_[0]  = 0;
      repxAccept_[1]  = 0;
      nParameters_ = system.perturbation().getNParameters();
   }

   /*
   * Destructor.
   */
   ReplicaMove::~ReplicaMove()
   {
      if (ptPositionPtr_) delete [] ptPositionPtr_;
      if (myPositionPtr_) delete [] myPositionPtr_;
   }

   /*
   * Read attempting interval.
   */
   void ReplicaMove::readParam(std::istream& in)
   {
      readBegin(in,"ReplicaMove");
      read<long>(in, "interval", interval_);
      if (interval_ <= 0) {
         UTIL_THROW("Invalid value input for interval_");
      }
      readEnd(in);
   }

   /*
   * Find the tempering variable and create data files.
   */
   void ReplicaMove::initialize()
   {
      // Get the perturbation parameter
      myParam_.allocate(nParameters_);
      ptParam_.allocate(nParameters_);
      for (int i = 0; i < nParameters_; ++i) {
         myParam_[i] = system().perturbation().parameter(i);
      }
      // Allocate memory
      int nAtom = system().simulation().atomCapacity();
      ptPositionPtr_ = new Vector[nAtom];
      myPositionPtr_ = new Vector[nAtom];
   }

   /*
   * Perform replica exchange move.
   */
   bool ReplicaMove::move()
   {
      MPI::Request request[8];
      MPI::Status  status;
      double myWeight, ptWeight, exponential;
      int    isLeft, iAccept, myPort, ptPort;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int iA;
      
      // Default value for no replica exchange
      isLeft  = -1;
      iAccept = 0;

      if ((myId_ < nProcs_ - 1) && (stepCount_ >= myId_) &&
         ((stepCount_ - myId_) % (nProcs_ - 1)  ==0)) {
         // we are exchanging with processor myId_ + 1
         isLeft = 1;
         ptId_ = myId_ + 1;

      } else if ((myId_ > 0) && (stepCount_ >= myId_ - 1) &&
                ((stepCount_  - (myId_ - 1)) % (nProcs_ -1 ) == 0)) {
         // we are exchanging with processor myId_ - 1
         isLeft = 0;
         ptId_ = myId_ - 1; 
      } 
      if (isLeft == 0 || isLeft == 1) {
         // Set the port value for message tag
         myPort = myId_%2;
         ptPort = ptId_%2;
 
         // Update accumulator
         repxAttempt_[isLeft] += 1;

         // Exchange coupling parameters with partner
         for (int i = 0; i < nParameters_; ++i) {
            request[0] = communicatorPtr_->Irecv(&ptParam_[i], 1, MPI::DOUBLE, ptId_,
                                              TagParam[ptPort]);
            request[1] = communicatorPtr_->Isend(&myParam_[i], 1, MPI::DOUBLE, ptId_,
                                              TagParam[myPort]);
         
            // Synchronizing
            request[0].Wait();
            request[1].Wait();
         }   
         myWeight = system().perturbation().difference(ptParam_);
         
         // Collect tempering weights and make decision
         if (isLeft == 1) {

            // Receive energy difference from the right box
            request[2] = communicatorPtr_->Irecv(&ptWeight, 1, MPI::DOUBLE,
                                                ptId_, TagDecision[ptPort]);
            request[2].Wait();

            exponential = exp(-myWeight - ptWeight);

         } else {

            // Send energy difference to the left box
            request[2] = communicatorPtr_->Isend(&myWeight, 1, MPI::DOUBLE,
                                                 ptId_, TagDecision[myPort]);
            request[2].Wait();
         }

         // Collect tempering weights and make decision
         if (isLeft == 1) {
   
            // Metropolis test
            iAccept = system().simulation().random().
                      metropolis(exponential) ? 1 : 0;

            // Send decision to the right box
            request[3] = communicatorPtr_->Isend(&iAccept, 1, MPI::INT,
                                                 ptId_, TagDecision[myPort]);
            request[3].Wait();

         } else {

            // Receive decision from the left box
            request[3] = communicatorPtr_->Irecv(&iAccept, 1, MPI::INT,
                                                 ptId_, TagDecision[ptPort]);
            request[3].Wait();

         }
         // Exchange particle configurations if the move is accepted
         if (iAccept == 1) {
         
            // Update accumulator
            repxAccept_[isLeft] += 1;

            Vector myBoundary;
            myBoundary = system().boundary().lengths();
                  
            Vector ptBoundary;
                  
            // Accomodate new boundary dimensions.
            request[4] = communicatorPtr_->Irecv(&ptBoundary, 1,
                                                 MpiTraits<Vector>::type, ptId_, TagConfig[ptPort]);

            // Send old boundary dimensions.
            request[5] = communicatorPtr_->Isend(&myBoundary, 1,
                                                  MpiTraits<Vector>::type, ptId_, TagConfig[myPort]);

            request[4].Wait();
            request[5].Wait();
                   
            system().boundary().setLengths(ptBoundary);

            // Pack atomic positions and types.
            iA = 0;
            for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
               for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
                  for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                     myPositionPtr_[iA] = atomIter->position();
                     iA++;
                  }
               }
            }
            
            // Accomodate new configuration.
            request[6] = communicatorPtr_->Irecv(ptPositionPtr_, iA,
                             MpiTraits<Vector>::type, ptId_, TagConfig[ptPort]);

            // Send old configuration.
            request[7] = communicatorPtr_->Isend(myPositionPtr_, iA,
                             MpiTraits<Vector>::type, ptId_, TagConfig[myPort]);

            request[6].Wait();
            request[7].Wait();
            
            // Adopt the new atomic positions.
            iA = 0;
            for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
               for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
                  for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter){
                     atomIter->position() = ptPositionPtr_[iA];
                     ++iA;
                  }
               }
            }
            

            // Notify component observers.
            //notifyObservers(ptId_);
            Notifier<int>::notifyObservers(ptId_);

         }  else {
            // the move was not accepted, do nothing
         }

      }


      // Output results needed to build exchange profile.
      if (isLeft == -1) {
        outputFile_ << "n" << std::endl;
      } else if ( isLeft != -1 && iAccept == 0) {
        outputFile_ << myId_ << std::endl;
      } else if ( isLeft != -1 && iAccept == 1) {
        outputFile_ << ptId_ << std::endl;
      }

      stepCount_++;

      return (iAccept == 1 ? true : false);

   }

}
#endif
#endif // ifdef UTIL_MPI
#endif // ifdef MCMD_PERTURB
