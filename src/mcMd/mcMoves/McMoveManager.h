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

      /// Constructor.
      McMoveManager(McSimulation& simulation);

      /// Destructor.
      ~McMoveManager();

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in parameter file input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Initialize at beginning of simulation run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      */
      McMove& chooseMove();

      /**
      * Return probability of move i.
      *
      * \param i index for McMove.
      * \return probability of move number i.
      */
      double probability(int i);

      /**
      * Output statistics for all moves.
      */
      void output();

      /**
      * Save state to a binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchive& ar);

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

   /// Return probability of move number i
   inline double McMoveManager::probability(int i)
   {
      assert(i >= 0);  
      assert(i < size());  
      return probabilities_[i];
   }

}
#endif
