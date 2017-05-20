#ifndef TOOLS_PAIR_ENERGY_H
#define TOOLS_PAIR_ENERGY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/analyzers/Analyzer.h>       // base class 
#include <tools/neighbor/CellList.h>        // member
#include <simp/interaction/pair/LJPair.h>   // member
#include <tools/neighbor/CellList.h>        // member
#include <util/containers/GArray.h>         // member

namespace Tools
{

   /**
   * Analyzer to compute total nonbonded pair energy.
   */
   class PairEnergy : public Analyzer
   {
   public:

      /**
      * Constructor.
      *
      * \param configuration  parent Configuration 
      */
      PairEnergy(Configuration &configuration);

      /**
      * Constructor.
      *
      * \param processor  parent Processor 
      */
      PairEnergy(Processor &processor);

      /**
      * Constructor.
      *
      * \param configuration  parent Configuration
      * \param fileMaster  associated FileMaster
      */
      PairEnergy(Configuration &configuration, FileMaster& fileMaster);

      /** 
      * Read parameters from file.
      *
      * \param in  input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Setup before main loop.
      */
      virtual void setup();

      /** 
      * Compute nonbonded pair energy.
      *
      * \param iStep  step counter
      */
      virtual void sample(long iStep);

      /**
      * Output final results to file.
      */
      virtual void output();

   private:

      /// Pair interaction type (hard-coded for now).
      typedef Simp::LJPair Interaction;

      // TODO: Generalize to other interaction types

      /// Output file stream
      std::ofstream outputFile_;
      
      /// Current interaction
      Interaction interaction_;

      /// CellList used to calculate energies
      CellList cellList_;

      /// Cutoff length for potential
      double cutoff_;

      /// Maximum number of atoms in cell list
      int atomCapacity_;

      /// Has readParam been called?
      int isInitialized_;

      /// Store timesteps this analyzer was run
      GArray<int> timesteps_;

      /// Store energies for the runs
      GArray<double> energies_;

   };

}
#endif
