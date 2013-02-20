#ifndef MCMD_MC_MOVE_CPP
#define MCMD_MC_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"
#include <mcMd/simulation/Simulation.h>
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McMove::McMove(Simulation& simulation) 
    : simulationPtr_(&simulation),
      randomPtr_(&simulation.random()),
      nAttempt_(0),
      nAccept_(0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   McMove::~McMove()
   {}

   /*
   * readParam, empty default implementation.
   */
   void McMove::readParameters(std::istream &in)
   {}
   
   /*
   * Read the probability from file.
   */
   void McMove::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }
   
   /*
   * Load internal state from an archive.
   */
   void McMove::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<double>(ar, "probability", probability_); 
      ar & nAttempt_;
      ar & nAccept_;
   }

   /*
   * Save internal state to an archive.
   */
   void McMove::save(Serializable::OArchive &ar)
   {  
      ar & probability_;
      ar & nAttempt_;
      ar & nAccept_;
   }

   /*
   * Trivial implementation - initializes counters.
   */
   void McMove::setup()
   { 
      nAttempt_ = 0;
      nAccept_  = 0;
   }

   /*
   * Trivial default implementation - always returns false.
   */
   bool McMove::move()
   { 
      ++nAttempt_;
      return false; 
   }

   /*
   * Trivial default implementation - do nothing
   */
   void McMove::output()
   {}

}
#endif
