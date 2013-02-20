#ifndef MCMD_MC_MOVE_H
#define MCMD_MC_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class Simulation;

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * The virtual move() method must generate a trial move, decide whether
   * to accept or reject it, and update the associated System or Systems
   * if it is accepted.
   *
   * \ingroup McMd_McMove_Module
   */
   class McMove : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      McMove(Simulation& simulation);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~McMove();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Set the probability for this McMove.
      *
      * \param probability Probability of this move being chosen. 
      */
      void setProbability(double probability);

      /**
      * Setup before the beginning of each simulation run.
      *
      * This method zeros the statistical accumulators nAttempt
      * and nAccept. Derived class implementations should 
      * complete any other initialization that is required at 
      * the beginning of each simulation run within a sequence
      * of several such runs.
      */
      virtual void setup();

      /**
      * Generate, attempt, and accept or reject a Monte Carlo move.
      *
      * Implementations of this by subclasses should:
      *
      *    - Generate a move (e.g., choose new atomic positions)
      *    - Decide to accept or reject (use random::metropolis)
      *    - Implement the move if accepted
      *
      * Trivial default implemention always returns false.
      *
      * \return true if accepted, false if rejected
      */
      virtual bool move();

      // Accessor Functions

      /**
      * Return probability for this McMove.
      */
      double probability() const;

      /**
      * Return number of moves that have been attempted.
      */
      long nAttempt() const;

      /**
      * Return number of moves that have been accepted.
      */
      long nAccept() const;

      /**
      * Output statistics for this move (called at the end of the simulation)
      */
      virtual void output();

   protected:

      /**
      * Increment the number of attempted moves.
      */
      void incrementNAttempt();

      /**
      * Increment the number of accepted moves.
      */
      void incrementNAccept();

      /**
      * Get parent Simulation object.
      */
      Simulation& simulation();

      /**
      * Get Random number generator of parent Simulation.
      */
      Random& random();

      /**
      * Read the probability from file.
      */
      void readProbability(std::istream& in);

   private:

      /// Pointer to parent Simulation object
      Simulation  *simulationPtr_;

      /// Pointer to random number generator
      Random  *randomPtr_;

      /// Probability of choosing this move
      double  probability_;

      /// Number of moves that have been attempted by this object.
      long  nAttempt_;

      /// Number of moves that have been accepted by this object.
      long  nAccept_;

   };

   // Public inline methods

   /*
   * Return number of moves that have been attempted.
   */
   inline long McMove::nAttempt() const
   {  return nAttempt_; }

   /*
   * Return number of moves that have been accepted.
   */
   inline long McMove::nAccept() const
   {  return nAccept_; }

   // Protected inline methods

   /*
   * Increment the number of attempted moves.
   */
   inline void McMove::incrementNAttempt()
   {  ++nAttempt_; }

   /*
   * Increment the number of accepted moves.
   */
   inline void McMove::incrementNAccept()
   {  ++nAccept_; }

   /*
   * Get parent Simulation object.
   */
   inline Simulation& McMove::simulation()
   {  return *simulationPtr_; }

   /*
   * Get Random number generator.
   */
   inline Random& McMove::random()
   {  return *randomPtr_; }

   /*
   * Get the probability.
   */
   inline double McMove::probability() const
   {  return probability_; }

   /*
   * Set the probability.
   */
   inline void McMove::setProbability(double probability)
   {  probability_ = probability; }

}
#endif
