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
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/random/Random.h>
#include <util/misc/Observer.h>

#include <sstream>
#include <cmath>

namespace McMd
{

   using namespace Util;

   // Define static constant MPI message Tag.

   /*
   * Constructor.
   */
   ReplicaMove::ReplicaMove(System& system) :
      systemPtr_(&system),
      communicatorPtr_(0),
      myId_(-1),
      nProcs_(0),
      outputFile_(),
      nParameters_(0),
      interval_(-1),
      nSampling_(-1),
      ptPositionPtr_(0),
      myPositionPtr_(0),
      swapAttempt_(0),
      swapAccept_(0)
   {
      // Precondition
      if (!system.hasPerturbation()) {
         UTIL_THROW("Parent System has no Perturbation");
      }

      setClassName("ReplicaMove");
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
   void ReplicaMove::readParameters(std::istream& in)
   {
      read<long>(in, "interval", interval_);
      if (interval_ <= 0) {
         UTIL_THROW("Invalid value input for interval_");
      }
      read<int>(in, "nSampling", nSampling_);
      if (nSampling_ <= 0) {
         UTIL_THROW("Invalid value input for nSampling_");
      }

      // Allocate memory
      int nAtom = system().simulation().atomCapacity();
      ptPositionPtr_ = new Vector[nAtom];
      myPositionPtr_ = new Vector[nAtom];
   }

   /*
   * Load internal state from an archive.
   */
   void ReplicaMove::loadParameters(Serializable::IArchive &ar)
   {
      // Load parameters
      loadParameter<long>(ar, "interval", interval_);
      loadParameter<int>(ar, "nSampling", nSampling_);
      ar & swapAttempt_;
      ar & swapAccept_;

      // Validate
      if (interval_ <= 0) {
         UTIL_THROW("Invalid value input for interval_");
      }
      if (nSampling_ <= 0) {
         UTIL_THROW("Invalid value input for nSampling_");
      }

      // Allocate memory
      int nAtom = system().simulation().atomCapacity();
      ptPositionPtr_ = new Vector[nAtom];
      myPositionPtr_ = new Vector[nAtom];
   }

   /*
   * Save internal state to an archive.
   */
   void ReplicaMove::save(Serializable::OArchive &ar)
   {
      ar & interval_;
      ar & nSampling_;
      ar & swapAttempt_;
      ar & swapAccept_;

      int nAtom = system().simulation().atomCapacity();
      ptPositionPtr_ = new Vector[nAtom];
      myPositionPtr_ = new Vector[nAtom];
   }
   
   /*
   * Perform replica exchange move.
   */
   bool ReplicaMove::move()
   {
      MPI::Request request[4];
      MPI::Status  status;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int iA;
      int recvPt, sendPt;
  
      DArray<int> permutation;
      permutation.allocate(nProcs_);

      // Gather all derivatives of the perturbation Hamiltonians and parameters on processor with rank 0
      DArray<double> myDerivatives;
      myDerivatives.allocate(nParameters_);
      DArray<double> myParameters;
      myParameters.allocate(nParameters_);

      for (int i=0; i< nParameters_; i++) {
         myDerivatives[i] = system().perturbation().derivative(i);
         myParameters[i] = system().perturbation().parameter(i);
      }

      int size = 0;
      size += memorySize(myDerivatives);
      size += memorySize(myParameters);

      if (myId_ != 0) {

         MemoryOArchive sendCurrent;
         sendCurrent.allocate(size);

         sendCurrent << myDerivatives;
         sendCurrent << myParameters;

         sendCurrent.send(*communicatorPtr_, 0);
      } else {
         DArray< DArray<double> > allDerivatives;
         DArray< DArray<double> > allParameters;
         allDerivatives.allocate(nProcs_);
         allParameters.allocate(nProcs_);
   
         allDerivatives[0].allocate(nParameters_);
         allDerivatives[0] = myDerivatives;
         allParameters[0].allocate(nParameters_);
         allParameters[0] = myParameters;

         for (int i = 1; i<nProcs_; i++) {
            MemoryIArchive recvPartner;
            recvPartner.allocate(size);
            recvPartner.recv(*communicatorPtr_, i);
            allDerivatives[i].allocate(nParameters_);
            allParameters[i].allocate(nParameters_);
            recvPartner >> allDerivatives[i];
            recvPartner >> allParameters[i];
         }

         // Now we have the complete matrix U_ij = u_i(x_j), permutate nsampling steps according
         // to acceptance criterium
  
         // start with identity permutation
         for (int i = 0; i < nProcs_; i++)
            permutation[i] = i;

         for (int n =0; n < nSampling_; n++) {
            swapAttempt_++;
            // choose a pair i,j, i!= j at random
            int i = system().simulation().random().uniformInt(0,nProcs_);
            int j = system().simulation().random().uniformInt(0,nProcs_-1);
            if (i<=j) j++;

            // apply acceptance criterium
            double weight = 0;
            for (int k = 0; k < nParameters_; k++) {
               double deltaDerivative = allDerivatives[i][k] - allDerivatives[j][k];
               // the permutations operate on the states (the perturbation parameters)
               weight += (allParameters[permutation[j]][k] - allParameters[permutation[i]][k])*deltaDerivative;
             }
            double exponential = exp(-weight);
            int accept = system().simulation().random(). metropolis(exponential) ? 1 : 0;
   
            if (accept) {
               swapAccept_++;
               // swap states of pair i,j
               int tmp = permutation[i];
               permutation[i] = permutation[j];
               permutation[j] = tmp;
               }
         }
    
         // send exchange partner information to all other processors
         for (int i = 0; i < nProcs_; i++) {
            if (i != 0)
                communicatorPtr_->Send(&permutation[i], 1, MPI::INT, i, 0);
            else
                sendPt = permutation[i];

            if (permutation[i] != 0)
               communicatorPtr_->Send(&i, 1, MPI::INT, permutation[i], 1);
            else
               recvPt = i;
         }
      }

      
      if (myId_ != 0) {
         // partner id to receive from
         communicatorPtr_->Recv(&sendPt, 1, MPI::INT, 0, 0);
         // partner id to send to
         communicatorPtr_->Recv(&recvPt, 1, MPI::INT, 0, 1);
      }

      if (recvPt == myId_ || sendPt == myId_) {
         // no exchange necessary
         outputFile_ << sendPt << std::endl;
         return true;
      }

      assert(recvPt != myId_ && sendPt != myId_);

      Vector myBoundary;
      myBoundary = system().boundary().lengths();
            
      Vector ptBoundary;
            
      // Accomodate new boundary dimensions.
      request[0] = communicatorPtr_->Irecv(&ptBoundary, 1,
                                           MpiTraits<Vector>::type, recvPt, 1);

      // Send old boundary dimensions.
      request[1] = communicatorPtr_->Isend(&myBoundary, 1,
                                            MpiTraits<Vector>::type, sendPt, 1);

      request[0].Wait();
      request[1].Wait();
             
      system().boundary().setOrthorhombic(ptBoundary);

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
      request[2] = communicatorPtr_->Irecv(ptPositionPtr_, iA,
                       MpiTraits<Vector>::type, recvPt, 2); 

      // Send old configuration.
      request[3] = communicatorPtr_->Isend(myPositionPtr_, iA,
                       MpiTraits<Vector>::type, sendPt, 2);

      request[2].Wait();
      request[3].Wait();
      
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
      sendRecvPair partners;
      partners[0] = sendPt;
      partners[1] = recvPt;
      Notifier<sendRecvPair>::notifyObservers(partners);

      // Log information about exchange partner to file
      outputFile_ << sendPt << std::endl;

      return true;

   }

}
#endif
#endif // ifdef UTIL_MPI
#endif // ifdef MCMD_PERTURB
