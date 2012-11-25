#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_REPLICA_MOVE_H 
#define MCMD_REPLICA_MOVE_H 


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids 
* 
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu).
* Gibbs sampler algorithm implemented by Jens Glaser, March 2012.
*
* Distributed under the terms of the GNU General Public License.  
*/ 

#include <util/param/ParamComposite.h> // base class 
#include <util/space/Vector.h>          // Util namespace
#include <util/misc/Notifier.h>          // Util namespace
#include <util/containers/DArray.h>
#include <util/containers/Pair.h>
#include <util/global.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   class System;

   /**
    * A pair of receiving and sending partner ranks
    */
   typedef Pair<int> sendRecvPair;


   /**
   * Replica exchange Monte Carlo move using Gibbs sampling
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
   * This move implements a version replica exchange that uses the Gibbs sampler
   * technique to sample permutations in replica space from the corresponding
   * distribution.  The sampling is done by running a small (inexpensive)
   * Markov chain Monte Carlo simulation of \b nSamling steps at every replica
   * exchange step, at the end of which a permutation is obtained that is used
   * to swap configurations between replicas. This implies that all
   * configurations are permuted simultaneously (but there may be
   * configurations which stay on the same processor).
   *
   * The technique is described in detail in
   * John D. Chodera and Michael R. Shirts, J. Chem. Phys. 135, 194110 (2011)
   * 
   * \ingroup McMd_Perturb_Module
   */
   class ReplicaMove : public ParamComposite,
                       public Notifier<sendRecvPair>
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
      * The parameter block takes the parameter \b interval, which 
      * sets the interval between successive replica exchange attempts.
      *
      * At every replica exchange attempt, the permutation is sampled
      * from a MCMC run of \b nSampling steps.
      * Empirically, \b nSampling should be on the order of P^3 .. P^5,
      * where P is the number of processors.
      *
      * \param in input stream from which to read parameters.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive& ar);

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
      * Notify observers of a successful replica exchange
      *
      * \param partners a pair of indices of partner replicas. Needs to be known
      *  for communication.
      */
      void notifyObservers(sendRecvPair partners);
 
      /**
      * Number of swap attempts
      */
      long nAttempt();
       
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

      /// Current processor's rank.
      int   myId_;

      /// Number of processors.
      int   nProcs_;

      /// Output file stream storing the acceptance statistics.
      std::ofstream outputFile_;

      /// Number of perturbation parameters.
      int   nParameters_;
      
      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Number of state swaps before exchanging
      int nSampling_;

      /// Pointer to allocated buffer to store atom positions.
      Vector   *ptPositionPtr_;

      /// Local copy of the system atom pointer.
      Vector   *myPositionPtr_;

      /// Count of attempted swaps
      long  swapAttempt_;

      /// Count of accepted swaps
      long  swapAccept_;

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
