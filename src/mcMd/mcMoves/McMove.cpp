#ifndef MC_MOVE_CPP
#define MC_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   void McMove::readParam(std::istream &in)
   {}
   
   /*
   * Read the probability from file.
   */
   void McMove::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }
   
   /*
   * Trivial implementation - initializes counters.
   */
   void McMove::initialize()
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
   {
   }

   /*
   * Save the internal state to an archive.
   */
   void McMove::save(Serializable::OArchiveType& ar)
   {  
      ar & nAttempt_;
      ar & nAccept_;
   }

   /**
   * Load the internal state to an archive.
   */
   void McMove::load(Serializable::IArchiveType& ar)
   {  
      ar & nAttempt_;
      ar & nAccept_;
   }

}
#endif
