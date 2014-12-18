#ifndef TOOLS_PAIR_ENERGY_AVERAGE_H
#define TOOLS_PAIR_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/analyzers/Analyzer.h>       // base class template
#include <tools/neighbor/CellList.h>        // member
#include <inter/pair/LJPair.h>              // member
#include <tools/neighbor/CellList.h>        // member
#include <util/containers/GArray.h>         // member

namespace Tools
{

   class PairEnergy : public Analyzer
   {

      // TODO: create Interaction class in tools for other interaction types
      typedef Inter::LJPair Interaction;

      public:
         /**
         * Constructor.
         *
         * \param configuration parent Configuration object
         */
         PairEnergy(Configuration &configuration);

         /**
         * Constructor.
         *
         * \param processor reference to parent Processor
         */
         PairEnergy(Processor &processor);

         /**
         * Constructor.
         *
         * \param configuration reference to parent Configuration
         * \param fileMaster reference to associated FileMaster
         */
         PairEnergy(Configuration &configuration, FileMaster& fileMaster);

         /** 
         * Read parameters from file.
         *
         * \param in input parameter stream
         */
         virtual void readParameters(std::istream& in);

         /** 
         * setup stuff
         */
         virtual void setup();

         /** 
         * Loop over cells and calculate energy
         *
         * \param iStep counter for number of steps
         */
         virtual void sample(long iStep);

         /**
         * Output results to file after simulation is completed.
         */
         virtual void output();

      private:
         /// Output file stream
         std::ofstream outputFile_;
         
         /**
         * Current Interaction
         */
         Interaction interaction_;

         /// CellList used to calculate energies
         CellList cellList_;

         /// Cutoff length for potential
         double cutoff_;

         /// atomCapacity for cell list
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
