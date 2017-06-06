#ifndef G1MSD_H
#define G1MSD_H

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/MeanSqDispArray.h>  // member template 
#include <util/space/Vector.h>                   // template parameter
#include <util/containers/DArray.h>             // member template

namespace Simp {
   class Species;
}
  
namespace McMd
{
  
   using namespace Util;  
   using namespace Simp;  

   /**
   * Autocorrelation for vector separation of any two monomers on a molecule.
   */
   class G1MSD : public SystemAnalyzer<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      G1MSD(System& system);
  
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
      DArray<Vector>    truePositions_;
   
      /// Array of position vectors, one per molecule of species.
      DArray<Vector>    oldPositions_;
   
      /// Array of shift vectors, one per molecule of species.
      DArray<IntVector> shifts_;
   
      /// Pointer to relevant Species.
      Species *speciesPtr_;
   
      /// Index of relevant Species.
      int     speciesId_;
   
      /// Number of molecules in the species (must not change).
      int     nMolecule_;
   
      /// Maximum length of each sequence in AutoCorrArray.
      int     capacity_;
      
      /// Number of atoms per molecule of this species.
      int     nAtom_;     
   
   };

}
#endif
