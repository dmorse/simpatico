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

      #if 0
      beginReadManager(in);

      initFactory();
 
      std::string name;
      McMove* mcMovePtr;

      // Loop over McMoves 
      bool isEnd = false;
      while (!isEnd) {

         readBlank(in);

         // Read and instantiate a new McMove
         mcMovePtr = factory().readObject(in, *this, name, isEnd);

         if (!isEnd) {
            if (mcMovePtr != 0) {
               append(*mcMovePtr, name);
            } else {
              UTIL_THROW("McMove subclass name not recognized"); 
            }
         }

      }


      // Add closing bracket to param output format
      End* endPtr = &addEnd();
      if (ParamComponent::echo() && isParamIoProcessor()) { 
         endPtr->writeParam(Log::file());
      }
      #endif

      Manager<McMove>::readParam(in);

      probabilities_.allocate(size());

      // Store and normalize probabilities
      double  totalProbability = 0.0;
      int     iMove;
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }

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

   /*
   * Save state to binary file archive.
   */
   void McMoveManager::save(Serializable::OArchive& ar)
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].save(ar);
      }
   }

   /*
   * Load state from a binary file archive.
   */
   void McMoveManager::load(Serializable::IArchive& ar)
   {
      for (int i=0; i < size(); ++i) {
         (*this)[i].load(ar);
      }
   }

}
#endif
