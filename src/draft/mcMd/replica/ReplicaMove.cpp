#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef REPLICA_MOVE_CPP
#define REPLICA_MOVE_CPP


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ReplicaMove.h"                      

#include <mcMd/perturb/Perturbation.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/misc/Observer.h>

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
      myParam_(0.0),
      systemPtr_(&system),
      communicatorPtr_(0),
      nProcs_(0),
      myId_(-1),
      ptId_(-1),
      xFlag_(0),
      repxAttempt_(),
      repxAccept_(),
      ptPositionPtr_(0),
      myPositionPtr_(0),
      outputFile_(),
      energyFile_(),
      interval_(-1)
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
      myParam_ = system().perturbation().parameter();

      // Allocate memory
      int nAtom = system().simulation().atomCapacity();
      ptPositionPtr_ = new Vector[nAtom];
      myPositionPtr_ = new Vector[nAtom];
   }

   /*
   * Empty implementation, to be filled by derived classes.
   */
   bool ReplicaMove::move()
   {
      MPI::Request request[8];
      MPI::Status  status;
      double ptParam,  myWeight, ptWeight;
      int    isLeft, iAccept, myPort, ptPort;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomPtr;
      int iA;

      // Default value for idle processor
      isLeft  = -1;
      iAccept = 0;

      // Idenfity active processor and its parter ID; left one has smaller ID.
      if ((myId_ + xFlag_)%2 == 0 && myId_ < nProcs_-1) {
         isLeft = 1;
         ptId_ = myId_ + 1;
      } else if ((myId_ + xFlag_)%2 == 1 && myId_ > 0) {
         isLeft = 0;
         ptId_ = myId_ - 1;
      }

      // Start to talk with partner
      if (isLeft == 1 || isLeft == 0) {

         // Set the port value for message tag
         myPort = myId_%2;
         ptPort = ptId_%2;

         // Update accumulator
         repxAttempt_[isLeft] += 1;

         // Exchange coupling parameters with partner
         request[0] = communicatorPtr_->Irecv(&ptParam, 1, MPI::DOUBLE, ptId_,
                                              TagParam[ptPort]);
         request[1] = communicatorPtr_->Isend(&myParam_, 1, MPI::DOUBLE, ptId_,
                                              TagParam[myPort]);

         // Synchronizing
         request[0].Wait();
         request[1].Wait();

         // Evaluating energy change.
         myWeight = system().perturbation().difference(ptParam);

         // Collect tempering weights and make decision
         if (isLeft == 1) {

            // Receive energy difference from the right box
            request[2] = communicatorPtr_->Irecv(&ptWeight, 1, MPI::DOUBLE,
                                                 ptId_, TagDecision[ptPort]);
            request[2].Wait();

            // Metropolis test
            iAccept = system().simulation().random().
                      metropolis(exp(-myWeight - ptWeight)) ? 1 : 0;

           // Output the two weights and the metropolis of their summation.
           // energyFile_ << setw(10) << myWeight << setw(10)
           //             << ptWeight << setw(2) << iAccept << std::endl;


            // Send decision to the right box
            request[3] = communicatorPtr_->Isend(&iAccept, 1, MPI::INT,
                                                 ptId_, TagDecision[myPort]);
            request[3].Wait();

         } else {

            // Send energy difference to the left box
            request[2] = communicatorPtr_->Isend(&myWeight, 1, MPI::DOUBLE,
                                                 ptId_, TagDecision[myPort]);
            request[2].Wait();

            // Receive decision from the left box
            request[3] = communicatorPtr_->Irecv(&iAccept, 1, MPI::INT,
                                                 ptId_, TagDecision[ptPort]);
            request[3].Wait();

         }

         // Exchange particle configurations if the move is accepted
         if (iAccept == 1) {
         
            // Update accumulator
            repxAccept_[isLeft] += 1;

            // Pack atomic positions and types.
            iA = 0;
            for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
               for (system().begin(iSpec, molIter); !molIter.isEnd(); ++molIter){
                  for (molIter->begin(atomPtr); !atomPtr.isEnd(); ++atomPtr) {
                     myPositionPtr_[iA] = atomPtr->position();
                     iA++;
                  }
               }
            }

            // Accomodate new configuration.
            request[4] = communicatorPtr_->Irecv(ptPositionPtr_, iA,
                             MpiTraits<Vector>::type, ptId_, TagConfig[ptPort]);

            // Send old configuration.
            request[5] = communicatorPtr_->Isend(myPositionPtr_, iA,
                             MpiTraits<Vector>::type, ptId_, TagConfig[myPort]);

            request[4].Wait();
            request[5].Wait();

            // Adopt the new atomic positions.
            iA = 0;
            for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
               for (system().begin(iSpec, molIter); !molIter.isEnd(); ++molIter){
                  for (molIter->begin(atomPtr); !atomPtr.isEnd(); ++atomPtr) {
                     atomPtr->position() = ptPositionPtr_[iA];
                     ++iA;
                  }
               }
            }

            // Notify component observers.
            //notifyObservers(ptId_);
            Notifier<int>::notifyObservers(ptId_);

         }  else {
         }

      }


      // Output results needed to build exchange profile.
      outputFile_ << ((isLeft != -1 && iAccept == 1) ? ptId_ : myId_)
                  << std::endl;

      // Flip the value of xFlag_ before exit.
      xFlag_ = (xFlag_ == 0 ? 1 : 0);

      return (iAccept == 1 ? true : false);

   }

}
#endif
#endif // ifdef UTIL_MPI
#endif // ifdef MCMD_PERTURB
