#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef REPLICA_MOVE_H 
#define REPLICA_MOVE_H 


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids 
* 
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu) 
* Distributed under the terms of the GNU General Public License.  
*/ 

#include <util/param/ParamComposite.h> // base class 
#include <util/space/Vector.h>          // Util namespace
#include <util/util/Notifier.h>          // Util namespace

#include <util/global.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Replica exchange Monte Carlo move.
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
   * In a replica exchange simulation, each replica lives on a different
   * MPI processor, with processors 0, 1, 2, ... ordered in order of
   * increasing values of the tempering parameter. In each attempted move,
   * each replica (possibly excluding one of the end replicas) attempts to 
   * exchange its configuration with one of its neighbors. The choice of 
   * which pairs of neighbors attempt to exchange is controlled by an int 
   * member xFlag_ whose value toggles between 0 and 1 after every attempt. 
   * When xFlag_ = 0, the following pairs attempt exchanges:
   *        0<->1, 2<->3, ...  (left replicas are ranked with even number)
   * When xFlag_ = 1, the other groups of pairs attempt exchanges:
   *        1<->2, 3<->4, ...  (left replicas are ranked with odd number)
   *
   * \ingroup Perturb_Module
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
      double myParam_;

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

      /// Flag determining which replica pairs attempt the exchange.
      int   xFlag_;

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

      /// Output file stream storing the acceptance statistics.
      std::ofstream energyFile_;

      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Tags for exchanging parameters.
      static const int TagParam[2];

      /// Tags for exchanging energy/decision.
      static const int TagDecision[2];

      /// Tags for exchanging configuration.
      static const int TagConfig[2];
      std::ofstream    repMoveFile_;

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
