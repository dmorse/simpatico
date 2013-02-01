#ifndef MCMD_MC_MOVE_MANAGER_CPP
#define MCMD_MC_MOVE_MANAGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <mcMd/mcMoves/McMoveManager.h>
#include <mcMd/mcMoves/McMoveFactory.h>
#include <mcMd/mcSimulation/McSimulation.h>

#include <util/random/Random.h>

namespace McMd
{

   using namespace Util;

   // Constructor.
   McMoveManager::McMoveManager(McSimulation& simulation)
   : Manager<McMove>(),
     simulationPtr_(&simulation),
     systemPtr_(&simulation.system()),
     randomPtr_(&simulation.random())
   {  setClassName("McMoveManager"); }

   // Destructor
   McMoveManager::~McMoveManager()
   {}

   /**
   * Return a pointer to a new McMoveFactory object.
   */
   Factory<McMove>* McMoveManager::newDefaultFactory() const
   {  return new McMoveFactory(*simulationPtr_, *systemPtr_); }

   /* 
   * Read instructions for creating objects from file.
   */
   void McMoveManager::readParam(std::istream &in)
   {
      Manager<McMove>::readParam(in);

      // Allocate and store probabilities
      probabilities_.allocate(size());
      double  totalProbability = 0.0;
      int     iMove;
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }

      // Allocate and store and normalize probabilities
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }
   }

   /*
   * Load internal state from an archive.
   */
   void McMoveManager::loadParameters(Serializable::IArchive &ar)
   {
      Manager<McMove>::loadParameters(ar);
      ar & probabilities_;
   }

   /*
   * Load internal state from an archive.
   */
   void McMoveManager::save(Serializable::OArchive &ar)
   {
      Manager<McMove>::save(ar);
      ar & probabilities_;
   }

   /*
   * Initialize all moves just prior to a run.
   */
   void McMoveManager::setup()
   {
      for (int iMove = 0; iMove < size(); ++iMove) {
         (*this)[iMove].setup();
      }
   }

   /*
   * Choose a McMove at random.
   */
   McMove& McMoveManager::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   void McMoveManager::output()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }

}
#endif
