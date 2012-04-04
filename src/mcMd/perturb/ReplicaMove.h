#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_REPLICA_MOVE_H 
#define MCMD_REPLICA_MOVE_H 


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids 
* 
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu) 
* Distributed under the terms of the GNU General Public License.  
*/ 

#include <util/param/ParamComposite.h> // base class 
#include <util/space/Vector.h>          // Util namespace
#include <util/util/Notifier.h>          // Util namespace
#include <util/containers/DArray.h>
#include <util/global.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Staggered Replica exchange Monte Carlo move.
   *
   * This class implements a rather general form of the algorithm for a 
   * replica exchange / parallel tempering move.  The algorithm assumes that 
   * the Boltzmann weight W(X,p) for a system in configuration X depends upon 
   * some floating point parameter p in a manner that is described by an
   * Perturbation object associated with the parent System.  The ReplicaMove 
   * constructor requires a parent System as a parameter, and requires that 
   * the System has an associated Perturbation object when the ReplicaMove
   * is constructed.
   *
   * This implements a cyclic replica exchange, where at the first call
   * a move between replicas 0<>1 is attempted, at the second call betweeen
   * replicas 1<>2,  ..., at the (N-1)th call between N-1 <> N, then 
   * 0<>1 again and so forth.
   * The replica exchange interval is set in the parameter
   * block (see the readParam() documentation) and determines
   * how many integration steps are performed between replica exchange
   * intervals. A full sweep over all replicas takes
   * (number of replicas -1) attempts.
   *
   * If the replica exchange for a given pair is successful, as determined
   * by the Metropolis criterium associated with the Perturbation, the complete
   * configurations of the pair are exchanged (not just their parameters).
   *
   * \ingroup McMd_Perturb_Module
   */
   class ReplicaMove : public ParamComposite,
                       public Notifier<int>
   {
   
   public:
   
      /**
      * Constructor. 
      */
      ReplicaMove(System& system);
   
      /**
      * Destructor. 
      */
      virtual ~ReplicaMove();

      /**
      * Read parameters.
      * 
      * The parameter block takes the parameter \b interval,
      * which sets the interval between successive replica exchange attempts.
      * A full sweep of all replicas is accomplished after
      * \b interval*(number of replicas - 1) integration steps.
      *
      * \param in input stream from which to read parameters.
      */
      virtual void readParam(std::istream& in);

      /**
      * Initialize variables needed by a replica move. Needs to be invoked
      * prior to any replica move.
      */
      virtual void initialize();

      /**
      * Attempt, and accept or reject a replica exchange move.
      */
      virtual bool move();

      /**
      * Get interval (number of steps between successives attempts).
      */
      long interval() const;
   
      /**
      * Return true iff counter is a multiple of the interval.
      *
      * \param counter simulation step counter
      */
      bool isAtInterval(long counter) const;
      
      /**
      * Number of swap attempts
      *
      * \param left index for direction of attempted exchange.
      */
      long nAttempt();
       
      void notifyObservers(int recvPt);

      /**
      * Number of accepted swaps
      */
      long nAccept(); 

   protected:

      /**
      * Return the associated system by reference.
      */
      System& system();

   private:

      /// System reference.
      System* systemPtr_;

      /// Get the communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;

      /// Number of processors.
      int   nProcs_;

      /// Current processor's rank.
      int   myId_;

      /// Number of perturbation parameters.
      int   nParameters_;
      
      /// Count of attempted swaps
      long  swapAttempt_;

      /// Count of accepted swaps
      long  swapAccept_;

      /// Number of state swaps before exchanging
      int nSampling_;

      /// Pointer to allocated buffer to store atom positions.
      Vector   *ptPositionPtr_;

      /// Local copy of the system atom pointer.
      Vector   *myPositionPtr_;

      /// Output file stream storing the acceptance statistics.
      std::ofstream outputFile_;

      /// Current number of steps
      long stepCount_;

      /// Number of simulation steps between subsequent actions.
      long interval_;

   };
   // Inline methods

   /*
   * Return interval value.
   */
   inline long ReplicaMove::interval() const
   { return interval_; }
   
   /*
   * Return true iff counter is a multiple of the interval.
   */
   inline bool ReplicaMove::isAtInterval(long counter) const
   { return (counter%interval_ == 0); }
   
   /*
   * Number of attempts in given direction.
   */
   inline long ReplicaMove::nAttempt() 
   {  return swapAttempt_; }
   
   /*
   * Number of accepted moves in given direction.
   */
   inline long ReplicaMove::nAccept()
   {  return swapAccept_; }

   /*
   * Return reference to parent System.
   */
   inline System& ReplicaMove::system()
   { 
      assert(systemPtr_);
      return *systemPtr_; 
   }

}
#endif
#endif // ifdef UTIL_MPI
#endif // ifdef MCMD_PERTURB
