#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_REPLICA_MOVE_H 
#define MCMD_REPLICA_MOVE_H 


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids 
* 
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.  
*/ 

#include <util/param/ParamComposite.h> // base class 
#include <util/space/Vector.h>          // Util namespace
#include <util/misc/Notifier.h>          // Util namespace
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
      * Number of attempts in specified direction.
      *
      * \param left index for direction of attempted exchange.
      */
      long nAttempt(int left);
       
      void notifyObservers(int partnerId);

      /**
      * Number of accepted moves in specified direction.
      *
      * \param left index for direction of attempted exchange.
      */
      long nAccept(int left); 

   protected:

      /**
      * Return the associated system by reference.
      */
      System& system();

   private:

      /// Tempering variable.
      DArray<double> myParam_;
      
      DArray<double> ptParam_;

      /// System reference.
      System* systemPtr_;

      /// Get the communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;

      /// Number of processors.
      int   nProcs_;

      /// Current processor's rank.
      int   myId_;

      /// Active neighboring (partner) replica's rank.
      int   ptId_;

      /// Number of perturbation parameters.
      int   nParameters_;
      
      /// Count the number of times the replica move is called to determine
      /// when this processor should attempt a replica exchange
      int   stepCount_;

      /// Count of attempted moves.
      long  repxAttempt_[2];

      /// Count of accepted moves.
      long  repxAccept_[2];

      /// Pointer to allocated buffer to store atom positions.
      Vector   *ptPositionPtr_;

      /// Local copy of the system atom pointer.
      Vector   *myPositionPtr_;

      /// Output file stream storing the acceptance statistics.
      std::ofstream outputFile_;


      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Tags for exchanging parameters.
      static const int TagParam[2];

      /// Tags for exchanging energy/decision.
      static const int TagDecision[2];

      /// Tags for exchanging configuration.
      static const int TagConfig[2];

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
   inline long ReplicaMove::nAttempt(int left) 
   {  return (left == 0 ? repxAttempt_[0]:repxAttempt_[1]); }
   
   /*
   * Number of accepted moves in given direction.
   */
   inline long ReplicaMove::nAccept(int left) 
   {  return (left == 0 ? repxAccept_[0]:repxAccept_[1]); }

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
