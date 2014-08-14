#ifndef SPAN_ATOM_MSD_H
#define SPAN_ATOM_MSD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/analyzers/Analyzer.h>            // base class template
#include <util/accumulators/MeanSqDispArray.h>  // member template 
#include <util/space/Vector.h>                  // member template parameter
#include <util/containers/DArray.h>             // member template

namespace SpAn
{

   using namespace Util;

   /**
   * Mean-squared displacement of specific atom in molecules of specific species.
   *
   * \ingroup SpAn_Analyzer_Module
   */
   class AtomMSD : public Analyzer
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param processor reference to parent Processor
      */
      AtomMSD(Processor &processor);
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
  
      /** 
      * Determine number of molecules and allocate memory.
      */
      virtual void setup();
   
      /** 
      * Evaluate end-to-end vectors of all chains, add to ensemble.
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

      /// Statistical accumulator
      MeanSqDispArray<Vector> accumulator_;

      /// Array of position vectors, one per molecule of species.
      DArray<Vector> truePositions_;
   
      /// Array of position vectors, one per molecule of species.
      DArray<Vector> oldPositions_;
   
      /// Array of shift vectors, one per molecule of species.
      DArray<IntVector> shifts_;
   
      /// Index of relevant Species.
      int speciesId_;
   
      /// Local index of atom1
      int atomId_;
   
      /// Number of molecules in the species (must not change).
      int nMolecule_;
   
      /// Maximum length of each sequence in AutoCorrArray.
      int capacity_;
   
      /// Has readParam been called?
      int isInitialized_;
   
   };

}
#endif
