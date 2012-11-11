#ifndef MCMD_MC_MOVE_MANAGER_H
#define MCMD_MC_MOVE_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Manager.h>          // base class template
#include "McMove.h"                      // base class template parameter
#include <util/containers/DArray.h>      // member template

namespace Util { class Random; }

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * Manager for a set of McMove objects.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_McMove_Module
   */
   class McMoveManager : public Manager<McMove>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent McSimulation
      */
      McMoveManager(McSimulation& simulation);

      /**
      * Destructor.
      */
      ~McMoveManager();

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream &in);

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
      * Initialize at beginning of simulation run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      *
      * \return chosen McMove
      */
      McMove& chooseMove();

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;

      /**
      * Output statistics for all moves.
      */
      void output();

   private:

      /// Array of McMove probabilities.
      DArray<double>  probabilities_;

      /// Pointer to random number generator.
      McSimulation* simulationPtr_;

      /// Pointer to random number generator.
      McSystem* systemPtr_;

      /// Pointer to random number generator.
      Random* randomPtr_;

      /// Return pointer to a new McMoveFactory.
      virtual Factory<McMove>* newDefaultFactory() const;

   };

   // Inline functions

   /*
   * Return probability of move number i
   */
   inline double McMoveManager::probability(int i) const
   {
      assert(i >= 0);  
      assert(i < size());  
      return probabilities_[i];
   }

}
#endif
